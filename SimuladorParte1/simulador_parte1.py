import numpy as np
import sys
import click

def fill_matrixs(Gn, In, content):
    for component in content:
        splitted_line = component.split(" ")
        
        if component[0] == "I":
            node_a = int(splitted_line[1])
            node_b = int(splitted_line[2])
            I = float(splitted_line[4])
            In[node_a] += -I
            In[node_b] += I

        elif component[0] == "R":
            node_a = int(splitted_line[1])
            node_b = int(splitted_line[2])
            R = float(splitted_line[3])
            Gn[node_a, node_a] += 1/R
            Gn[node_b, node_b] += 1/R
            Gn[node_a, node_b] += -1/R
            Gn[node_b, node_a] += -1/R
            
        elif component[0] == "G":
            node_a = int(splitted_line[1])
            node_b = int(splitted_line[2])
            node_c = int(splitted_line[3])
            node_d = int(splitted_line[4])
            G = float(splitted_line[5])
            Gn[node_a, node_c] += G
            Gn[node_a, node_d] += -G
            Gn[node_b, node_c] += -G
            Gn[node_b, node_d] += G
    return Gn, In

def get_number_of_nodes(content):
    # Store only unique elements
    nodes = set()
    
    content_lenght = len(content)
    i = 0
    while(i < content_lenght):
        # Skip and removes lines with:       
        if content[i][0] == "\n" or content[i][0] == "*" or content[i][0] == " "\
        or content[i][0] == "\0": 
            
            content.remove(content[i])

            # The next element becomes the current one, so doesn't need to add 1 in "i".
            # Content lenght down by 1 
            content_lenght -= 1

        elif content[i][0] == "I" or content[i][0] == "R" or content[i][0] == "G":
            splitted_line = content[i].split(" ")
            
            # Add node 1 and 2
            nodes.add(int(splitted_line[1]))
            nodes.add(int(splitted_line[2]))

            if content[i][0] == "G":
                # Add control nodes 3 and 4
                nodes.add(int(splitted_line[3]))
                nodes.add(int(splitted_line[4]))
            i += 1
        else:
            sys.exit(f"Netlist com a {i}ª linha inválida")

    return len(nodes), content

def print_results(e):    
    print()
    for i in range(len(e)):
        print(f"e{i + 1} = {e[i]} V")
    print()

@click.command()
@click.option("--netlist", help="Netlist file path", prompt=True)
def main(netlist):
    """
    Calculates the nodal voltages of the circuit provided by the netlist.
    """
    
    # Get the content of the netlist file
    try:
        with open(netlist) as netfile:
            content = netfile.readlines()
    except:
        sys.exit("Error. Invalid netlist path.")

    # Get number of nodes for creating Gn and In
    n, content = get_number_of_nodes(content) 
    Gn = np.zeros((n, n))
    In = np.zeros(n)

    # Fills Gn and In according with the netlist 
    Gn, In = fill_matrixs(Gn, In, content)

    # Removes node 0 (GND)
    Gn = Gn[1:, 1:]
    In = In[1:]

    # Solve the system Gn.e = In and get "e"
    e = np.linalg.solve(Gn, In)
    
    print_results(e)

main()




