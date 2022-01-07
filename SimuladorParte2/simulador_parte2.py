import numpy as np
import sys
import click
from math import pi
from cmath import rect

w = 0

def Add_Idc(splitted_line, In):
    node_a = int(splitted_line[1])
    node_b = int(splitted_line[2])
    I = float(splitted_line[4])
    In[node_a] += -I
    In[node_b] += I
    return In

def Add_Isin(splitted_line, In):
    node_a = int(splitted_line[1])
    node_b = int(splitted_line[2])
    
    # Take care if Isin doesn't has a DC value space 
    if len(splitted_line) == 8:
        I_module = float(splitted_line[5])
        I_phase = float(splitted_line[7])*pi/180 # pass to radians
    
    else: # So, len(splitted_line) == 7
        I_module = float(splitted_line[4])
        I_phase = float(splitted_line[6])*pi/180 # pass to radians
    
    In[node_a] += -rect(I_module, I_phase)
    In[node_b] += rect(I_module, I_phase)
    return In

def Add_R(splitted_line, Yn):
    node_a = int(splitted_line[1])
    node_b = int(splitted_line[2])
    R = float(splitted_line[3])
    
    Yn[node_a, node_a] += 1/R
    Yn[node_b, node_b] += 1/R
    Yn[node_a, node_b] += -1/R
    Yn[node_b, node_a] += -1/R
    return Yn

def Add_GmIdc(splitted_line, Yn):
    node_a = int(splitted_line[1])
    node_b = int(splitted_line[2])
    node_c = int(splitted_line[3])
    node_d = int(splitted_line[4])
    G = float(splitted_line[5])
    
    Yn[node_a, node_c] += G
    Yn[node_a, node_d] += -G
    Yn[node_b, node_c] += -G
    Yn[node_b, node_d] += G
    return Yn

def Add_C(splitted_line, Yn):
    node_a = int(splitted_line[1])
    node_b = int(splitted_line[2])
    C = float(splitted_line[3])
    
    Yn[node_a, node_a] += w*C*1j
    Yn[node_b, node_b] += w*C*1j
    Yn[node_a, node_b] += -w*C*1j
    Yn[node_b, node_a] += -w*C*1j
    return Yn

def Add_L(splitted_line, Yn):
    node_a = int(splitted_line[1])
    node_b = int(splitted_line[2])
    L = float(splitted_line[3])
    
    Yn[node_a, node_a] += 1/(w*L*1j)
    Yn[node_b, node_b] += 1/(w*L*1j)
    Yn[node_a, node_b] += -1/(w*L*1j)
    Yn[node_b, node_a] += -1/(w*L*1j)
    return Yn

def Add_K(splitted_line, Yn):
    node_a = int(splitted_line[1])# = 0
    node_b = int(splitted_line[2])# = 1
    node_c = int(splitted_line[3])# = 0
    node_d = int(splitted_line[4])# = 2
    L1 = float(splitted_line[5])
    L2 = float(splitted_line[6])
    M = float(splitted_line[7])

    L = np.array([[L1, M],
                  [M, L2]])
    r = np.linalg.inv(L)
    estampa = np.array([[ r[0][0]/(w*1j),  -r[0][0]/(w*1j),   r[0][1]/(w*1j),  -r[0][1]/(w*1j)],
                        [-r[0][0]/(w*1j),   r[0][0]/(w*1j),  -r[0][1]/(w*1j),   r[0][1]/(w*1j)],
                        [ r[1][0]/(w*1j),  -r[1][0]/(w*1j),   r[1][1]/(w*1j),  -r[1][1]/(w*1j)],
                        [-r[1][0]/(w*1j),   r[1][0]/(w*1j),  -r[1][1]/(w*1j),   r[1][1]/(w*1j)]], dtype="complex_")
    
    print(estampa)

    nodes = [node_a, node_b, node_c, node_d]
    for node_1 in nodes:
        for node_2 in nodes:
            Yn[node_1, node_2] += estampa[nodes.index(node_1), nodes.index(node_2)]

    return Yn

def fill_matrixs(Yn, In, content):
    for component in content:
        splitted_line = component.split(" ")
        
        if component[0] == "I":
            if splitted_line[3] == "DC":
                In = Add_Idc(splitted_line, In)
            elif splitted_line[3] == "SIN":
                In = Add_Isin(splitted_line, In)

        elif component[0] == "R":
            Yn = Add_R(splitted_line, Yn)
            
        elif component[0] == "G":
            Yn = Add_GmIdc(splitted_line, Yn)

        elif component[0] == "C":
            Yn = Add_C(splitted_line, Yn)
        
        elif component[0] == "L":
            Yn = Add_L(splitted_line, Yn)

        elif component[0] == "K":
            Yn = Add_K(splitted_line, Yn)

    return Yn, In

def get_number_of_nodes(content):
    """
    Get number of nodes for creating the matrices

    Args:
        content (list): Full content of netlist file in "readlines" format.

    Returns:
        [int]: Number of nodes.
        [list]: Content of netlist file with useless line removed.
    """
    
    # Store only unique elements
    nodes = set()
    
    i = 0
    content_lenght = len(content)
    while(i < content_lenght):
        # Skip and removes blank lines and comments:       
        if content[i][0] == "\n" or content[i][0] == "*" or content[i][0] == " "\
           or content[i][0] == "\0": 
            
            content.remove(content[i])

            # The next element becomes the current one, so doesn't need to add 1 in "i".
            # Content lenght down by 1 
            content_lenght -= 1

        elif content[i][0] == "I" or content[i][0] == "R" or content[i][0] == "G"\
            or content[i][0] == "L" or content[i][0] == "C" or content[i][0] == "K":
            
            splitted_line = content[i].split(" ")
            
            # Add node "a" and "b"
            nodes.add(int(splitted_line[1]))
            nodes.add(int(splitted_line[2]))

            if content[i][0] == "G" or content[i][0] == "K":
                # nodes "c" and "d" (control nodes for "G", and nodes for the secondary 
                # inductor of the transformer "K")
                nodes.add(int(splitted_line[3]))
                nodes.add(int(splitted_line[4]))

            # Get frequency of sine wave source
            elif splitted_line[3] == "SIN":
                global w
                w = 2*pi*float(splitted_line[6])

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

    # Get number of nodes for creating Yn and In
    n, content = get_number_of_nodes(content) 
    Yn = np.zeros((n, n), dtype = 'complex_')
    In = np.zeros(n, dtype = 'complex_')

    # Fills Yn and In according with the netlist 
    Yn, In = fill_matrixs(Yn, In, content)

    # Removes node 0 (YnD)
    Yn = Yn[1:, 1:]
    In = In[1:]

    # Solve the system Yn.e = In and get "e"
    e = np.linalg.solve(Yn, In)
    
    print_results(e)

main()