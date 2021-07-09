"""

This script allows for the creation of networks that visualizes similarity
of different PDB structures.

"""
import requests     # To make REST API calls to Protein Data Bank API.
import json         # To interact with JSON response from Protein Data Bank.

import pandas as pd
import networkx as nx
from pyvis.network import Network

# Helper function for query_structure_similarity_pdbs()
def query_structure_similarity_pdbs_json_processor(json_dic):

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

# PDB search function that gets PDB ID and structure similarity score for pdb_id_to_query
def query_structure_similarity_pdbs(pdb_id_to_query, mode='strict_shape_match', return_all=True):

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
    
    Return
    ------
    Returns a list of tuples, each containing a PDB ID and structure similarity score
    of PDB structures that are similar to the PDB whose ID was queried.

    If PDB ID queried does not exist, the function will return None
    
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
    if (response.ok) and (response.status_code == 200):
        return query_structure_similarity_pdbs_json_processor(response.json())
    else:
        return None

def load_to_dataframe(similar_pdbs):
    """
    Returns a Panda DataFrame containing node and edge data.

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

# Uses a PDB ID as a "seed" to get structure similarity data
def get_structure_similarity_data_by_seed(seed_id, mode='strict_shape_match', max_search_depth=1):
    """
    Using a seed PDB ID, an in-depth structure similarity search is conducted up until the
    limit set by max_search_depth.

    Parameters
    ----------
    seed_id: str
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
            search_id = seed_id
            
            similar_structures = query_structure_similarity_pdbs(search_id, mode=mode)
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
                    
                    similar_structures = query_structure_similarity_pdbs(search_id, mode=mode)
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


# Helper function for get_structure_similarity_data_from_range()
# Modified method from pdbGenerator.py in PR#1
def get_pdb_ids_in_range(lower, upper=None):
    """
    Returns a list of strings representing PDB IDs in the range of [lower, upper]
    If there are issues with the inputs, the function will return False

    Parameters
    ----------
    lower : str
        A PDB ID that will be the lower bound of our search.
        The lower range input will be considered a "wildcard" if the string contains
        an '*' after a value (e.g. 1A*) or the PDB ID is incomplete (e.g. 1A).
    
    upper : str
        A PDB ID that will be the upper bound of our search.
        The upper range input will be considered a "wildcard" if the string contains
        an '*' after a value (e.g. 1A*) or the PDB ID is incomplete (e.g. 1A).
        If upper == None, then range is from [lower, 9ZZZ].
    
    Corner Cases
    ------------
    Input Errors
    1. lower or upper are not of type str
    2. lower > upper; that is, lower comes after upper (e.g. BBBB to AAAA).
    3. The first character in the string for either lower or upper is a letter; not proper PDB ID format.
    4. len(lower) or len(upper) is greater than 4; not proper PDB ID format
    5. lower and upper contain upper case letters; just convert lower and upper to lower case.
    6. Keep in Mind: wildcard symbol or a space is found in between two sets of letters (e.g. 1A*B, *BBB, " AA" or "AD C")
    
    Special Ranges
    1. lower == upper; in this case, the function will attempt to download 1 file.
    2. lower == '*' and upper == '*'; in this case, we will download all files.
    3. lower is a valid ID, but upper == '*'; in this case, we will download all files between lower and upper == '9ZZZ'
    4. lower == '*', but upper is a valid ID; in this case, we will download all files between lower == '0000' and upper
    
    """
        
    alpha = ['a','b','c','d','e','f','g','h','i','j','k','l','m',
                'n','o','p','q','r','s','t','u','v','w','x','y','z']
    num = ['0','1','2','3','4','5','6','7','8','9']
    
    # Idea behind having the numbers go before the letters
    # is to count up from 0 to Z. So 0000 would be the first
    # ID, 0001 would be the second, and so on
    numalpha = num + alpha


    ### Checking for Default Values; specifically if upper == None or dir == None.
    if upper == None:
        # Set upper to '9zzz', so that range is [lower, '9zzz']
        upper = '9zzz'
    else:
        pass

    

    ### Check to see if lower, upper, and dir are string
    
    # Checking lower and upper
    if (isinstance(lower, str) == False) or (isinstance(upper, str) == False):    # Input Error 1: lower and upper are not str.
        return False
    else:
        # Input Error 5: lower and upper contain upper case letters; just convert lower and upper to lower case.
        lower = lower.lower()
        upper = upper.lower()


    ### lower, upper, and dir are confirmed strings.


    ### Input Error 4: len(lower) or len(upper) is greater than 4; not proper PDB ID format
    if (len(lower) > 4) or (len(upper) > 4):                  
        return False
    else:
        pass


    ### Time to preprocess them for wildcards before continuing with more input error checks.

    # Checking for Special Ranges
    # If lower == upper, our loop code should be able to handle that case on its own.
    if lower == '*' and upper == '*': # Special Range 2: lower == '*' and upper == '*', so we will download all files
        lower = '1000'
        upper = '9zzz'
    elif upper == '*':                # Special Range 3: lower is a valid ID, but upper == '*'
        upper = '9zzz'
    elif lower == '*':                # Special Range 4: lower == '*', but upper is a valid ID
        lower = '1000'
    else:
        pass

    # Preprocess inputs for wildcards
    if '*' in lower:
        # If lower contains a '*', that means a wildcard is present.
        # Find the '*', and fill in the remaining characters with '0'
        for idx in range(len(lower)):
            if lower[idx] == '*':
                lower = lower[0: idx] + (numalpha[0] * (4 - idx))
                break
            else:
                pass

    if len(lower) < 4:
        # If len(lower) < 4, ID contains a wildcard through incompletion.
        # Because length is less than 4, we will add the remaining characters as '0'.
        lower = lower + (numalpha[0] * (4 - len(lower)))

    if '*' in upper:
        # If upper contains a '*', that means a wildcard is present.
        # Find the '*', and fill in the remaining characters with '0'
        for idx in range(len(upper)):
            if upper[idx] == '*':
                upper = upper[0: idx] + (numalpha[0] * (4 - idx))
                break
            else:
                pass
    
    if len(upper) < 4:
        # If len(upper) < 4, ID contains a wildcard through incompletion.
        # Because length is less than 4, we will add the remaining characters as '0'.
        upper = upper + (numalpha[0] * (4 - len(upper)))


    ### Continue with Input Error checks
    if (lower > upper) and (lower != '*') and (upper != '*'):   # Input Error 2: lower > upper (MAY BE SOURCE OF BUGS)
        return False
    elif (lower[0] in alpha) or (upper[0] in alpha):            # Input Error 3: The first character in either lower or upper is a letter
        return False
    elif (lower[0] == num[0]) or (upper[0] == num[0]):          # Input Error: The first character is '0'; not proper PDB format
        return False
    else:
        pass
    

    ### Inputs have been preprocessed. Time to get the PDB files
    lower_firstpos = numalpha.index(lower[0])
    lower_secondpos = numalpha.index(lower[1])
    lower_thirdpos = numalpha.index(lower[2])
    lower_fourthpos = numalpha.index(lower[3])

    currentID = ''

    ###### List to return
    created_ids = []

    for firstpos in range (lower_firstpos, 10):
        
        if currentID > upper:
            # currentID is greater than upper, stop loop
            break   
        
        for secondpos in range(lower_secondpos, 36):
            
            if currentID > upper:
                # currentID is greater than upper, stop loop
                break   
            
            for thirdpos in range(lower_thirdpos, 36):
                
                if currentID > upper:
                    # currentID is greater than upper, stop loop
                    break   
                
                for fourthpos in range(lower_fourthpos, 36):
                    
                    currentID = numalpha[firstpos] + numalpha[secondpos] + numalpha[thirdpos] + numalpha[fourthpos]

                    if currentID > upper:
                        # currentID is greater than upper, stop loop
                        break  
                    else:
                        created_ids.append(currentID)


    # Return the list of IDs
    return created_ids

# Uses a range of PDBs to get structure similarity data
def get_structure_similarity_data_from_range(lower, upper=None, mode='strict_shape_match', num_neighbors=10):
    
    """
    Searches for structure similarity for PDBs whose IDs are in the range [lower, upper]

    lower: str
        Lower end PDB ID for range.

    upper: str
        Upper end PDB ID for range.

    mode: str
        Structure similarity search mode to use. Either 'strict_shape_match' or 'relaxed_shape_match'

    num_neighbors: int
        Maximum number of similar structures to include in graph for each ID.

    Return
    ------
    A list of tuples with the following structure: (PDB ID searched, list of similar PDBS of searched ID)
    If num_neighbors is not an int or it is less than 0, it will return None.
    
    """
    # Input check
    # lower and upper will have their "input check" within structure_similarity_search()

    if (isinstance(num_neighbors, int) == False) or (num_neighbors < 0):
        return None
    elif isinstance(mode, str) == False:
        return None
    else:
        pass

    # Step 1: Get all possible PDB IDs in range
    possible_pdb_ids = get_pdb_ids_in_range(lower=lower, upper=upper)

    # Step 2: Get structure similarity data for each ID
    results = []

    for id in possible_pdb_ids:

        # For debug
        print('Querying: ' + id)

        data = query_structure_similarity_pdbs(id, mode=mode)

        # Check that we got data
        if data == None:
            
            # For debug
            print(id + ' contains no data')
            
            continue
        else:
            
            # For debug
            print(id + ' contains data')


        # Trim data array based on num_neighbors
        # and create the tuple to append to "results" list

        data.pop(0)     # Remember to pop the first element as it is the ID that was queried

        if len(data) >= num_neighbors:
            data = data[:10]
        else:
            pass

        tup = (id, data)
        results.append(tup)

    return results


def main():

    # Parameters
    lower = '1A00'
    upper = '1A10'
    search_mode = 'strict_shape_match'      # Mode used in structure similarity search; either 'strict_shape_match' or 'relaxed_shape_match'
    num_neighbors = 10

    # Step 1: Retrieve structure similarity data used for graph.
    similar_pdbs = get_structure_similarity_data_from_range(lower=lower, upper=upper, mode=search_mode, num_neighbors=num_neighbors)
    
    # Step 2: Load the data into a pandas data frame to make it easier to create nodes and edges.
    df = load_to_dataframe(similar_pdbs)

    # Step 3: Load the DataFrame data into a NetworkX graph.
    graph = nx.from_pandas_edgelist(df, source='Source', target='Target', edge_attr='Weight')

    # Step 4: Load the NetworkX graph into pyvis
    net = Network(height='1000px', width='1000px')
    net.from_nx(graph)

    physics_params = {
        'gravity': 0,
        'central_gravity': 0,
        'spring_length': 0,
        'spring_strength': 0,
        'damping': 0,
        'overlap': 10
    }

    net.options.physics.use_barnes_hut(physics_params)
    
    # Step 5: Display network
    net.show('Pyvis_test.html')

    # For debug
    print()
    print('Number of nodes: ' + str(len(net.nodes)))
    print('Number of edges: ' + str(len(net.edges)))
   

if __name__ == "__main__":
    main()