"""

This script creates a Dash app that displays a visdcc network that
displays the relationship of structure similarity between more than one
PDB structure.

"""
import requests     # To make REST API calls to Protein Data Bank API.
import json         # To interact with JSON response from Protein Data Bank.
import uuid         # To make unique IDs for edges in case there are repeated edges.

import dash
import visdcc
import pandas as pd
import dash_html_components as html

import timeit         # For timing data processing step.

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

    If similar_pdbs is empty or None, function will return None

    Parameters
    ----------
        similar_pdbs: list of tuples with the following structure -> (PDB ID searched: str, similar PDBS for searched ID: list(tuple(ID, Similarity Score)))
            A list of tuples, with the above representation.

    """
    # If similar_pdbs is empty or None, return None
    if similar_pdbs == None:
        return None
    elif len(similar_pdbs) == 0:
        return None
    else:
        pass
    
    # Initialize DataFrame
    df = pd.DataFrame(data=None, index=None, columns=['Source', 'Target', 'Weight'], copy=None)

    for tup in similar_pdbs:
        source_id = tup[0]         # The PDB ID that was use in the search is consider the source
        
        # The list of PDBs associated with the search are consider the targets
        list_of_pdbs = tup[1]
        for target in list_of_pdbs:
            target_id = target[0]
            target_sim_score = target[1]

            new_row = {'Source': source_id, 'Target': target_id, 'Weight': target_sim_score}
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
            'id': str(uuid.uuid4()),
            'from': source,
            'to': target,
            'width': 1
        }

        edges.append(edge)

    return edges



def structure_similarity_search(starter_id, mode='strict_shape_match', max_search_depth=1):
    """
    Using a "starter" PDB ID, an in-depth structure similarity search is conducted up until the
    limit set by max_search_depth.

    Parameters
    ----------
    starter_id: str
        A PDB ID that is used to start the search.

    mode: str
        Mode to use when searching for similar structures.
        Valid modes: 'strict_shape_match', 'relaxed_shape_match'

    max_search_depth: int
        The maximum number of "searches" that can be done. One search means iterating
        through an entire list of PDBs.
    
    Return
    ------
    A list of tuples with the following structure: (PDB ID searched, list of similar PDBS of searched ID)
    If max_search_depth is not an int or it is less than 0, it will return None.
    """

    # Input Checks, specifically for max_num and max_search_depth as they will be used in our loops
    if isinstance(max_search_depth, int) == False:
        return None
    elif max_search_depth < 0:
        return None
    else:
        pass


    # Results matrix is a 2-D array that will be used for deep searches.
    # n x 1 matrix, so each row contains one item.
    results_matrix = []

    # Loop to do the in-depth search.
    for i in range(0, max_search_depth):
        
        if (len(results_matrix) == 0):
            # Do an initial search
            search_id = starter_id
            
            similar_structures = get_similar_pdbs(search_id, mode=mode)
            similar_structures.pop(0)   # First item in the list is the pdb represented by "id = current_id"
            
            tup = (search_id, similar_structures)
            
            results_matrix.append([tup])
        else:
            list_to_go_through = results_matrix[i - 1]
            list_to_add = []

            # In this for-loop, we use the fact that each list in the matrix is a tuple of (PDB ID, Similar Structures)
            for item in list_to_go_through:
                for tup in item[1]:
                    search_id = tup[0]
                    
                    similar_structures = get_similar_pdbs(search_id, mode=mode)
                    similar_structures.pop(0)   # First item in the list is the pdb represented by "id = current_id"

                    tup = (search_id, similar_structures)

                    list_to_add.append(tup)

            results_matrix.append(list_to_add)

    
    # Final result matrix.
    results = []

    for row in results_matrix:
        for col in results_matrix:
            for tup in col:
                results.append(tup)
    
    return results


def main():

    # Parameters
    pdbID = '5S5Y'                          # PDB ID to search for similar PDBs
    search_mode = 'strict_shape_match'      # Mode used in structure similarity search; either 'strict_shape_match' or 'relaxed_shape_match'
    search_depth = 2                        # How many layers of "search" to conduct.

    # Step 1: Create the Dash app that will be used to display our network
    app = dash.Dash()

    ###### Start our timer
    start_time = timeit.default_timer()
    print('Started timing')


    # Step 2: Retrieve structure similarity data used for graph.
    similar_pdbs = structure_similarity_search(pdbID, mode=search_mode, max_search_depth=search_depth)
    
    # Step 3: Load the data into a pandas data frame to make it easier to create nodes and edges.
    df = load_to_dataframe(similar_pdbs)
    
    # Step 4: Create nodes.
    nodes = create_nodes(df)

    ##### For debug, print number of nodes in graph
    print('Number of nodes = ' + str(len(nodes)))
    
    # Step 5: Create edges.
    edges = create_edges(df)



    ###### Stop our timer and print out time elapsed.
    stop_time = timeit.default_timer()
    print('Time Elapsed in Data Processing: ' + (str(stop_time - start_time)))


    
    # Step 6: Create visdcc Network and add it to a Dash layout.
    if nodes != None and edges != None:
        network = visdcc.Network(id = 'net',
                                data = {'nodes': nodes, 'edges': edges},
                                options = {
                                    'height': '600px',
                                    'width': '100%',
                                    'physics': {
                                        'enabled': False
                                    }
                                })

        app.layout = html.Div([network])
    else:
        app.layout = html.Div([html.H1('Error In Creating Graph')])
    
    # Step 7: Define any Dash callbacks.


    # Step 8: Start the app.
    app.run_server(debug=False)


if __name__ == "__main__":
    main()
