"""

Example of visualizing a small network of
10 PDB files, which are structurally similar,
using visdcc and Dash.

"""
import requests     # To make REST API calls to Protein Data Bank API
import json         # To interact with JSON response from Protein Data Bank

import dash
import visdcc
import pandas as pd
import dash_html_components as html

def get_similar_pdbs_process_json(json_dic):

    """
    A helper function for get_similar_pdbs().
    Processes the raw JSON from an HTTP response by
    extracting the PDB IDs and similarity scores, putting
    them in tuples, and appending them to a list.

    Parameters
    ----------
    json_dic: dictionary
        A dictionary that represents a JSON.

    Things to get from JSON:
        1. "result_set" - Array containing the search results.
            -> In this array, we want "identifier" to get PDB ID and "original_score" to get similarity score.
    """
    
    # Get "result_set" array
    # This is an array of dictionaries for each PDB
    results = json_dic['result_set']

    # Array filled with tuples of (PDB ID, Structure Similarity Score)
    ids_and_scores = []

    # Fill ids_and_scores with the tuples (PDB ID, Structure Similarity Score)
    for pdb in results:
        ID = pdb['identifier']
        score = pdb['services'][0]['nodes'][0]['original_score']

        ids_and_scores.append((ID, score))

    return ids_and_scores

def get_similar_pdbs(pdb_id_to_query, mode='strict_shape_match', return_all=True):

    """
    Returns a list of tuples, each containing a PDB ID and structure similarity score
    of PDB structures that are similar to the PDB whose ID was queried.

    Tuples in list are structured as: (PDB ID, Structure Similarity Score)

    Parameters
    ----------
    pdb_id_to_query: str
        PDB ID of structure to find similar structures.
    
    mode: str
        Mode to use when searching for similar structures.
        Valid modes: 'strict_shape_match', 'relaxed_shape_match'

    return_all: bool
        If true, returns all results from Protein Data Bank.
        If false, returns a few of the top results from Protein Data Bank.
    
    Input Checks
    ------------
    Naturally, there would be some code to check to see if the inputs of the
    function are valid. In this case, input checks may be unnecessary as bad
    inputs would result in an Error Code from the HTTP response.
    """

    # REST API endpoint
    url = 'https://search.rcsb.org/rcsbsearch/v1/query?json='

    # Dictionary representation of the JSON query
    query = {
        "query": {
            "type": "terminal",
            "service": "structure",
            "parameters": {
	            "operator": mode,
                "value": {
                    "entry_id": pdb_id_to_query,
                    "assembly_id": "1"
                }
            }
        },
        "request_options": {
            "return_all_hits": return_all
        },
        "return_type": "entry"
    }

    # Convert dictionary to a JSON string
    query_json = json.dumps(query)

    # Make the response request
    response = requests.get(url=(url + query_json))

    # If we get a code 200 (OK), return the array of tuples.
    if (response.ok):
        return get_similar_pdbs_process_json(response.json())
    else:
        return None

def load_to_dataframe(similar_pdbs):
    """
    Returns a Panda DataFrame that will be used to create nodes and edges.

    If similar_pdbs is empty, it will return None

    Parameters
    ----------
        similar_pdbs: list(tuple(PDB ID, Structure Similarity Score))
            A list of tuples, with representation (PDB ID, Structure Similarity Score), which
            represent a PDB structure that is structurally similar to the PDB at similar_pdbs[0]

    """
    # If similar_pdbs is empty, return None
    if len(similar_pdbs) == 0:
        return None
    
    # Initialize DataFrame
    df = pd.DataFrame(data=None, index=None, columns=['Source', 'Target', 'Weight'], copy=None)

    # First tuple in similar_pdbs contains the PDB ID of the structure that we were
    # looking for similar structures for.
    source = similar_pdbs[0]

    # Creates our nodes and our edge pairs.
    # Edges go from the 'Source' node (which is similar_pdbs[0]) to the 'Target' node.
    # The edges are weighted, with values equal to the similarity score.
    for target in similar_pdbs:
        
        if (target == source):
            # We don't want to have a circular edge.
            continue
        else:
            source_id = source[0]
            target_id = target[0]
            similarity_score = target[1]

            new_row = {'Source': source_id, 'Target': target_id, 'Weight': similarity_score}
            df = df.append(new_row, ignore_index=True)

    return df

def create_nodes(df):

    """
    Returns a list of dictionaries containing information about the nodes in the graph.

    Parameters
    ----------
    df : pandas.DataFrame
        Contains data placed into two columns: 'Source' and 'Target.'
        If df == None, we will return None
    
    """
    if df is None:
        return None

    # Makes a list of nodes from df
    node_list = list(set(df['Source'].unique().tolist() + df['Target'].unique().tolist()))
    
    # Creates the Nodes dictionary using node_list.
    # 'id' is the ID of the node.
    # 'label' is the label of the node.
    # 'shape' is how the node looks like when plotted.
    # 'size' is how big the node is.
    nodes = [{'id': node_name, 'label': node_name, 'shape': 'dot', 'size': 7}
                for i, node_name in enumerate(node_list)]

    return nodes

def create_edges(df):

    """
    Returns a list of dictionaries containing information about the weighted edges in the graph.

    Parameters
    ----------
    df : pandas.DataFrame
        Contains data placed into three columns: 'Source', 'Target', and 'Weight'.
        If df == None, we will return None

    """
    if df is None:
        return None

    edges = []

    for row in df.to_dict(orient='records'):
        source, target, score = row['Source'], row['Target'], row['Weight']
        
        edge = {
            'id': source + '__' + target,
            'from': source,
            'to': target,
            'width': 1
        }

        edges.append(edge)

    return edges


def main():

    # Example Constants
    pdbID = '5S5Y' # PDB ID to search for similar PDBs

    # Step 1: Create the Dash app that will be used to display our network
    app = dash.Dash()

    # Step 2: Retrieve 10 (or more) PDB files that are structurally similar to "pdbID".
    # The function returns a list of tuples
    similar_pdbs = get_similar_pdbs(pdbID)
    
    # Step 3: Load the data into a pandas data frame to make it easier to create nodes and edges.
    df = load_to_dataframe(similar_pdbs)
    
    # Step 4: Create nodes.
    nodes = create_nodes(df)
    
    # Step 5: Create edges.
    edges = create_edges(df)
    
    # Step 6: Create visdcc Network and add it to a Dash layout.
    if nodes != None and edges != None:
        network = visdcc.Network(id = 'net',
                                data = {'nodes': nodes, 'edges': edges},
                                options = {
                                    'height': '600px',
                                    'width': '100%'
                                })

        app.layout = html.Div([network])
    else:
        app.layout = html.Div([html.H1('Error In Creating Graph')])
    
    # Step 7: Define any Dash callbacks.


    # Step 8: Start the app.
    app.run_server(debug=False)


if __name__ == "__main__":
    main()
