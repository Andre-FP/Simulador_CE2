import numpy as np
import sys
import click
from math import pi
from cmath import rect

class Circuit:

    def __init__(self):
        self.n = None
        self.Gn = None
        self.In = None
        self.e = None
        self.w = 0
        self.components = []
        self.Idc = []
        self.Isin = []
        self.R = []
        self.GmI = []
        self.C = []
        self.L = []
        self.K = []
        self.vars = []
        self.vars_values = []
        self.i_values = []

    def newGn(self, Gn, newline=1):
        """Adds lines and columns to Gn matrix

        Args:
            Gn (np.ndarray): Conductance matrix
            newline (int, optional): Number of new lines. Defaults to 1.

        Returns:
            np.ndarray: new Gn matrix
        """ 
        n = len(Gn)
        newGn = np.zeros((n + newline, n + newline))
        newGn[:n, :n] = Gn[:, :]
        return newGn

    def newIn(self, In, newline=1):
        """Adds lines to In matrix

        Args:
            In (np.ndarray): Current vector
            newline (int, optional): Number of new lines. Defaults to 1.

        Returns:
            np.ndarray: new In matrix
        """
        n = len(In)
        newIn = np.zeros((n + newline, 1))
        newIn[:n] = In[:]
        return newIn

    def Add_VDC(self, splitted_line, Gn, In):
        """Add DC Voltage source to matrix equation 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """

        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        
        Gn = self.newGn(Gn)
        In = self.newIn(In)
        V = float(splitted_line[4])
        
        # Stores that a variable of interest was added to equation
        self.vars.append(f"j{splitted_line[0]}")

        In[-1] += -V

        Gn[-1, node_a] += -1         
        Gn[-1, node_b] +=  1
        Gn[node_a, -1] +=  1         
        Gn[node_b, -1] += -1

        return Gn, In

    def Add_E(self, splitted_line, Gn, In):
        """Add Voltage source coltrolled by voltage 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        node_c = int(splitted_line[3])
        node_d = int(splitted_line[4])

        Gn = self.newGn(Gn)
        In = self.newIn(In)
        
        # Stores that a variable of interest was added to equation
        self.vars.append(f"jx{splitted_line[0]}")

        Av = float(splitted_line[5])
        
        Gn[-1, node_a] +=  -1         
        Gn[-1, node_b] +=   1
        Gn[-1, node_c] +=  Av         
        Gn[-1, node_d] += -Av
        Gn[node_a, -1] +=   1         
        Gn[node_b, -1] +=  -1

        return Gn, In

    def Add_F(self, splitted_line, Gn, In):
        """Add Current source controlled by current 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        node_c = int(splitted_line[3])
        node_d = int(splitted_line[4])

        Gn = self.newGn(Gn)
        In = self.newIn(In)
        
        B = float(splitted_line[5])

        # Stores that a variable of interest was added to equation
        self.vars.append(f"jx{splitted_line[0]}")
        
        Gn[-1, node_c] += -1         
        Gn[-1, node_d] +=  1
        Gn[node_a, -1] +=  B         
        Gn[node_b, -1] += -B
        Gn[node_c, -1] +=  1         
        Gn[node_d, -1] += -1

        return Gn, In

    def Add_H(self, splitted_line, Gn, In):
        """Add Voltage source controlled by current 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        node_c = int(splitted_line[3])
        node_d = int(splitted_line[4])

        Gn = self.newGn(Gn, newline=2)
        In = self.newIn(In, newline=2)
        
        # Stores that two variables of interest was added to equation
        self.vars.append(f"jx{splitted_line[0]}")
        self.vars.append(f"jy{splitted_line[0]}")
        
        Rm = float(splitted_line[5])
        
        Gn[-2, node_c] += -1         
        Gn[-2, node_d] +=  1
        Gn[-1, node_a] += -1         
        Gn[-1, node_b] +=  1
        Gn[node_a, -1] +=  1         
        Gn[node_b, -1] += -1
        Gn[node_c, -2] +=  1         
        Gn[node_d, -2] += -1
        Gn[-1, -2]     += Rm          

        return Gn, In

    def Add_Idc(self, splitted_line, In):
        """Add DC Current source to matrix equation 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        I = float(splitted_line[4])
        In[node_a] += -I
        In[node_b] +=  I
        
        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "DC Current Source",
                     "Value": f"{I} A",
                     "Node a": node_a,
                     "Node b": node_b}
        self.Idc.append(meta_info)
        self.components.append(meta_info)

        return In

    def Add_Isin(self, splitted_line, In):
        """Add DC Voltage source to matrix equation 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        
         
        if len(splitted_line) == 8:
            I_module = float(splitted_line[5])
            I_phase = float(splitted_line[7])*pi/180 # pass to radians
        
        # Take care if Isin doesn't have a DC value space
        elif len(splitted_line) == 7: # So, len(splitted_line) == 7
            I_module = float(splitted_line[4])
            I_phase = float(splitted_line[6])*pi/180 # pass to radians
        
        else:
            sys.exit(f'\nInvalid netlist. Sine Current Source "{splitted_line[0]}" has an invalid number of arguments')

        In[node_a] += -rect(I_module, I_phase)
        In[node_b] += rect(I_module, I_phase)

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Sine Current Source",
                     "Module": f"{I_module} A",
                     "Phase": f"{I_phase*180/pi} grades",
                     "Node a": node_a,
                     "Node b": node_b}
        self.Isin.append(meta_info)
        self.components.append(meta_info)
        return In

    def Add_R(self, splitted_line, Gn):
        """Add Resistor to matrix equation 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        R = float(splitted_line[3])
        
        Gn[node_a, node_a] += 1/R
        Gn[node_b, node_b] += 1/R
        Gn[node_a, node_b] += -1/R
        Gn[node_b, node_a] += -1/R

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Resistor",
                     "Value": f"{R} Ω",
                     "Node a": node_a,
                     "Node b": node_b}
        self.R.append(meta_info)
        self.components.append(meta_info)
        return Gn

    def Add_GmI(self, splitted_line, Gn):
        """Add Current source controlled by voltage 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        node_c = int(splitted_line[3])
        node_d = int(splitted_line[4])
        G = float(splitted_line[5])
        
        Gn[node_a, node_c] += G
        Gn[node_a, node_d] += -G
        Gn[node_b, node_c] += -G
        Gn[node_b, node_d] += G

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Current Source Controlled by voltage",
                     "Value": f"{G} 1/A",
                     "Node a": node_a,
                     "Node b": node_b,
                     "Node c": node_c,
                     "Node d": node_d}
        self.GmI.append(meta_info)
        self.components.append(meta_info)
        return Gn

    def Add_C(self, splitted_line, Gn):
        """Add Capacitor to matrix equation 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        C = float(splitted_line[3])
        
        Gn[node_a, node_a] += self.w*C*1j
        Gn[node_b, node_b] += self.w*C*1j
        Gn[node_a, node_b] += -self.w*C*1j
        Gn[node_b, node_a] += -self.w*C*1j

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Capacitor",
                     "Value": f"{C} F",
                     "Node a": node_a,
                     "Node b": node_b}
        self.C.append(meta_info)
        self.components.append(meta_info)
        return Gn

    def Add_L(self, splitted_line, Gn):
        """Add Inductor to matrix equation 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        L = float(splitted_line[3])
        
        Gn[node_a, node_a] += 1/(self.w*L*1j)
        Gn[node_b, node_b] += 1/(self.w*L*1j)
        Gn[node_a, node_b] += -1/(self.w*L*1j)
        Gn[node_b, node_a] += -1/(self.w*L*1j)
        
        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Inductor",
                     "Value": f"{L} H",
                     "Node a": node_a,
                     "Node b": node_b}
        self.L.append(meta_info)
        self.components.append(meta_info)
        return Gn

    def Add_K(self, splitted_line, Gn):
        """Add Transformer to matrix equation 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
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
        estampa = np.array([[ r[0][0]/(self.w*1j),  -r[0][0]/(self.w*1j),   r[0][1]/(self.w*1j),  -r[0][1]/(self.w*1j)],
                            [-r[0][0]/(self.w*1j),   r[0][0]/(self.w*1j),  -r[0][1]/(self.w*1j),   r[0][1]/(self.w*1j)],
                            [ r[1][0]/(self.w*1j),  -r[1][0]/(self.w*1j),   r[1][1]/(self.w*1j),  -r[1][1]/(self.w*1j)],
                            [-r[1][0]/(self.w*1j),   r[1][0]/(self.w*1j),  -r[1][1]/(self.w*1j),   r[1][1]/(self.w*1j)]], dtype="complex_")
        
        nodes = [node_a, node_b, node_c, node_d]
        for node_1 in nodes:
            for node_2 in nodes:
                Gn[node_1, node_2] += estampa[nodes.index(node_1), nodes.index(node_2)]

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Transformer",
                     "Primary inductance": f"{L1} H",
                     "Secondary inductance": f"{L2} H",
                     "Mutual inductance": f"{M} H",
                     "Node a": node_a,
                     "Node b": node_b,
                     "Node c": node_c,
                     "Node d": node_d}
        self.L.append(meta_info)
        self.components.append(meta_info)
        return Gn


    def get_number_of_nodes(self, content):
        """
        Get number of nodes for creating the Gn and In matrices.
        Eliminates comments and blank lines.

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
                or content[i][0] == "L" or content[i][0] == "C" or content[i][0] == "K"\
                or content[i][0] == "F" or content[i][0] == "E" or content[i][0] == "H"\
                or content[i][0] == "V":

                splitted_line = content[i].split(" ")
                
                # Add node "a" and "b"
                nodes.add(int(splitted_line[1]))
                nodes.add(int(splitted_line[2]))

                if content[i][0] == "G" or content[i][0] == "K" or content[i][0] == "F"\
                    or content[i][0] == "E" or content[i][0] == "H":

                    # nodes "c" and "d" (control nodes for "G", and nodes for the secondary 
                    # inductor of the transformer "K")
                    nodes.add(int(splitted_line[3]))
                    nodes.add(int(splitted_line[4]))

                # Get frequency of sine wave source
                elif splitted_line[3] == "SIN":
                    self.w = 2*pi*float(splitted_line[6])        
                i += 1
            
            else:
                sys.exit(f"\nNetlist com a {i}ª linha inválida")

        #if self.w == 0:
        #    sys.exit("\nThe circuit must have at least one Sine Power Source. There aren't any in this netlist.")

        self.n = len(nodes)
        return self.n, content

    def fill_matrixs(self, Gn, In, content):
        """Calls each "Add" function to add each component to the
        KCL equation. 

        Args:
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector
            content (list): Full content of netlist file in "readlines" format.

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        for component in content:
            splitted_line = component.split(" ")
            
            if component[0] == "I":
                if splitted_line[3] == "DC":
                    In = self.Add_Idc(splitted_line, In)
                elif splitted_line[3] == "SIN":
                    In = self.Add_Isin(splitted_line, In)

            elif component[0] == "R":
                Gn = self.Add_R(splitted_line, Gn)
                
            elif component[0] == "G":
                Gn = self.Add_GmI(splitted_line, Gn)

            elif component[0] == "C":
                Gn = self.Add_C(splitted_line, Gn)
            
            elif component[0] == "L":
                Gn = self.Add_L(splitted_line, Gn)

            elif component[0] == "K":
                Gn = self.Add_K(splitted_line, Gn)
            
            elif component[0] == "F":
                Gn, In = self.Add_F(splitted_line, Gn, In)
            
            elif component[0] == "E":
                Gn, In = self.Add_E(splitted_line, Gn, In)

            elif component[0] == "H":
                Gn, In = self.Add_H(splitted_line, Gn, In)
            
            elif component[0] == "V":
                Gn, In = self.Add_VDC(splitted_line, Gn, In)

        return Gn, In

    def analyze(self, netlist):
        """Calculates the nodal voltages and variables of interest. 
        Stores values in attributes "e", "i_values", "vars_values".

        Args:
            netlist (str): Path for netlist file.

        Returns:
            np.ndarray: nodal voltages vector 
        """
        # Refresh all variables
        self.__init__()

        # Get the content of the netlist file
        try:
            with open(netlist) as netfile:
                content = netfile.readlines()
        except:
            sys.exit("Error. Invalid netlist path.")

        # Get the number of nodes for creating Gn and In
        n, content = self.get_number_of_nodes(content) 
        Gn = np.zeros((n, n))
        In = np.zeros((n, 1))

        # Fills Gn and In accordingly to the netlist 
        Gn, In = self.fill_matrixs(Gn, In, content)

        # Removes node 0 (GND)
        Gn = Gn[1:, 1:]
        In = In[1:]
        self.Gn = Gn
        self.In = In

        # Solve the system Gn.e = In and get "e"
        vars_values = np.linalg.solve(Gn, In)
        
        # Pick how many nodal voltage variables the system has  
        len_e = len(vars_values) - len(self.vars)

        # Saves the "e" values in self.e
        self.e = vars_values[:len_e] 

        # Saves i_values in self.i_values
        self.i_values = vars_values[len_e:]

        # Saves all variables in self.vars_values
        self.vars_values = vars_values
        return self.e

    def show(self, var):    
        """Show the value(s) of the variable of interest in screen.

        Args:
            var (str): Variable of interest.
                Possible values:
                - "Gn"
                - "In"
                - "vars"
                - "n"
                - "w"
                - "components"
                - "Idc"
                - "Isin"
                - "R"
                - "GmI"
                - "C"
                - "L"
                - "K"
        """
        var_stored = ""
        print("\n\n------------------------------------------------------------")

        if var == "Gn":
            print("\tAdmittance Matrix (Gn)")
            print(self.Gn) 
            
        elif var == "In":
            print("\tCurrent Nodes Matrix (In)")
            print(self.In) 
            
        elif var == "vars":    
            print("\tNodal Voltages and Currents")
            for i in range(len(self.e)):
                print(f"e{i + 1} = {self.e[i][0]:.3f} V")
            
            for i in range(len(self.i_values)):
                print(f"{self.vars[i]} = {self.i_values[i][0]:.3f} A")

            print("------------------------------------------------------------")
        
        elif var == "n":
            print(f"Number of nodes: {self.n}")

        elif var == "w":
            print(f"Circuit frequency: {self.w} rad/s")

        elif var == "components":
            var_stored = self.components
            title = "Circuit Components"
        elif var == "Idc":
            var_stored = self.Idc
            title = "DC Current Sources"
        elif var == "Isin":
            var_stored = self.Isin
            title = "Sines Current Sources"
        elif var == "R":
            var_stored = self.R
            title = "Circuit Resistors"
        elif var == "GmI":
            var_stored = self.GmI
            title = "Current Sources Controlled by Voltages"
        elif var == "C":
            var_stored = self.C
            title = "Circuit Capacitors"
        elif var == "L":
            var_stored = self.L
            title = "Circuit Inductors"
        elif var == "K":
            var_stored = self.K
            title = "Circuit Transformers"

        if var_stored != "":
            print(f"\t{title}")
            for info in var_stored:
                for atr, value in info.items():
                    print(f"{atr}: {value:.3f}")
                print()
            
@click.command()
@click.option("--netlist", help="Netlist file path", prompt=True)
def main(netlist):
    """
    Calculates the nodal voltages of the circuit provided by the netlist.
    """
    
    ckt = Circuit()
    ckt.analyze(netlist)
    ckt.show("vars")

main()