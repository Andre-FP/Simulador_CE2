from cmath import inf
import matplotlib.pyplot as plt
import numpy as np
import sys
import click
from math import pi, cos, exp
import random

class Circuit:

    def __init__(self):
        self.n = None
        self.Gn = None
        self.In = None
        self.e = []
        self.e0 = [] # Initial condition of nodal voltages (time t - deltaT)
        self.w = 0
        self.components = []
        self.Idc = []
        self.Isin = []
        self.Vsin = []
        self.R = []
        self.GmI = []
        self.C = []
        self.L = []
        self.K = []
        self.vars = []        # Name variables of interests (currents)
        self.vars_values = []
        self.i_values = []
        self.i_values0 = [] # Initial condition of i variables (time t - deltaT)
        self.i_vars0 = {} # Initial condition of i variables associating name and values
        self.deltaT = None 
        self.period = None

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

    def Add_Vsin(self, splitted_line, Gn, In, time):
        """Add AC Voltage source to matrix equation 

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
            V_dc = float(splitted_line[4])
            V_amplitude = float(splitted_line[5])
            V_w = float(splitted_line[6])*2*pi # pass to radians
            V_phase = float(splitted_line[7])*pi/180 # pass to radians

        # Take care if Vsin doesn't have a DC value space
        elif len(splitted_line) == 7: # So, len(splitted_line) == 7
            V_dc = 0
            V_amplitude = float(splitted_line[4])
            V_w = float(splitted_line[5])*2*pi # pass to radians
            V_phase = float(splitted_line[6])*pi/180 # pass to radians

        else:
            sys.exit(f'\nInvalid netlist. Sine Voltage Source "{splitted_line[0]}" has an invalid number of arguments')

        # Add voltage source current as variable of matrix equation (add 1 line to Gn and In)
        Gn, In = self.add_lines_matrixs(Gn, In, splitted_line[0], newline=1)

        Gn[-1, node_a] += -1         
        Gn[-1, node_b] +=  1
        Gn[node_a, -1] +=  1         
        Gn[node_b, -1] += -1

        In[-1] += - V_dc - V_amplitude*cos(V_w*time + V_phase)

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Sine Voltage Source",
                     "Module": f"{V_amplitude} A",
                     "Phase": f"{V_phase*180/pi} grades",
                     "DC value": f"{V_dc} A",
                     "Node a": node_a,
                     "Node b": node_b}
        self.Vsin.append(meta_info)
        self.components.append(meta_info)
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

    def Add_Isin(self, splitted_line, In, time):
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
            I_dc = float(splitted_line[4])
            I_amplitude = float(splitted_line[5])
            I_w = float(splitted_line[6])*2*pi # pass to radians
            I_phase = float(splitted_line[7])*pi/180 # pass to radians

        # Take care if Isin doesn't have a DC value space
        elif len(splitted_line) == 7: # So, len(splitted_line) == 7
            I_dc = 0
            I_amplitude = float(splitted_line[4])
            I_w = float(splitted_line[5])*2*pi # pass to radians
            I_phase = float(splitted_line[6])*pi/180 # pass to radians

        else:
            sys.exit(f'\nInvalid netlist. Sine Current Source "{splitted_line[0]}" has an invalid number of arguments')

        In[node_a] += - I_dc - I_amplitude*cos(I_w*time + I_phase)
        In[node_b] += I_dc + I_amplitude*cos(I_w*time + I_phase)

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Sine Current Source",
                     "Module": f"{I_amplitude} A",
                     "Phase": f"{I_phase*180/pi} grades",
                     "DC value": f"{I_dc} A",
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

    def Add_C(self, splitted_line, Gn, In):
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
        
        # Add capacitor current as variable of matrix equation (add 1 line to Gn and In)
        Gn, In = self.add_lines_matrixs(Gn, In, splitted_line[0], newline=1)

        # Resistor (Linearized capacitor)
        R = self.deltaT/(2*C)

        Gn[node_a, -1] +=  1
        Gn[node_b, -1] += -1
        Gn[-1, node_a] += -1
        Gn[-1, node_b] +=  1
        Gn[-1, -1]     +=  R

        # Get v(t-deltaT) and i(t-deltaT)
        if len(self.e0) == 0:
            v0 = float(splitted_line[4])
            i0 = 0
        else:
            if node_a == 0:
                eA = 0
            else:
                eA = self.e0[node_a - 1]
            
            if node_b == 0:
                eB = 0
            else:
                eB = self.e0[node_b - 1]

            v0 = eA - eB 
            i0 = self.i_vars0[f"j{splitted_line[0]}"]
            
        # Current source (Linearized capacitor)
        In[-1] += -v0 -R*i0 

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Capacitor",
                     "Value": f"{C} F",
                     "Node a": node_a,
                     "Node b": node_b}
        self.C.append(meta_info)
        self.components.append(meta_info)
        return Gn, In

    def add_lines_matrixs(self, Gn, In, component_name_id, newline=1):
        Gn = self.newGn(Gn, newline=newline)
        In = self.newIn(In, newline=newline)

        # Stores that a variable of interest was added to equation
        self.vars.append(f"j{component_name_id}")
        return Gn, In

    def Add_L(self, splitted_line, Gn, In):
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
        
        # Resistor (Linearized inductor)
        R = 2*L/(self.deltaT)

        # Add inductor current as variable of matrix equation (add 1 line to Gn and In)
        Gn, In = self.add_lines_matrixs(Gn, In, splitted_line[0], newline=1)

        Gn[node_a, -1] +=  1
        Gn[node_b, -1] += -1
        Gn[-1, node_a] += -1
        Gn[-1, node_b] +=  1
        Gn[-1, -1]     +=  R

        # Get v(t-deltaT) and i(t-deltaT)
        if len(self.e0) == 0:
            v0 = 0
            i0 = float(splitted_line[4])
        else:
            if node_a == 0:
                eA = 0
            else:
                eA = self.e0[node_a - 1]
            
            if node_b == 0:
                eB = 0
            else:
                eB = self.e0[node_b - 1]

            v0 = eA - eB 
            i0 = self.i_vars0[f"j{splitted_line[0]}"]

        # Current source (Linearized inductor)
        In[-1] += R*i0 + v0

        # Registering
        name = splitted_line[0]
        meta_info = {"Name": name,
                     "Component": "Inductor",
                     "Value": f"{L} H",
                     "Node a": node_a,
                     "Node b": node_b}
        self.L.append(meta_info)
        self.components.append(meta_info)

        return Gn, In

    def Add_K(self, splitted_line, Gn, In):
        """Add Transformer to matrix equation 

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
        L1 = float(splitted_line[5])
        L2 = float(splitted_line[6])
        M = float(splitted_line[7])

        #
        Gn = self.newGn(Gn, newline=2)
        In = self.newIn(In, newline=2)

        # Stores that two variables of interest was added to equation
        iab_name = f"j{splitted_line[0]}_{node_a}{node_b}"
        icd_name = f"j{splitted_line[0]}_{node_c}{node_d}"
        self.vars.append(iab_name)
        self.vars.append(icd_name)

        Gn[node_a, -2] +=  1
        Gn[node_b, -2] += -1
        Gn[node_c, -1] +=  1
        Gn[node_d, -1] += -1
        Gn[-2, node_a] += -1
        Gn[-2, node_b] += -1
        Gn[-2, -2]     += 2*L1/self.deltaT
        Gn[-2, -1]     += 2*M/self.deltaT
        Gn[-1, -2]     += 2*M/self.deltaT
        Gn[-1, -1]     += 2*L2/self.deltaT

        # Get v(t-deltaT) and i(t-deltaT)
        if len(self.e0) == 0:
            vab0 = 0
            vcd0 = 0
            iab0 = 0
            icd0 = 0
        else:
            nodes = [node_a, node_b, node_c, node_d]
            e = []
            for node in nodes:
                if node == 0:
                    e.append(0)
                else:
                    e.append(self.e0[node - 1])
            
            vab0 = e[0] - e[1]
            vcd0 = e[2] - e[3]
            iab0 = self.i_vars0[iab_name]
            icd0 = self.i_vars0[icd_name]

        In[-2] += (2/self.deltaT)*(L1*iab0 + M*icd0) + vab0
        In[-1] += (2/self.deltaT)*(M*iab0 + L2*icd0) + vcd0

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
        return Gn, In

    def Add_Diode(self, splitted_line, Gn, In, guess_v0diode):
        """Add Diode to matrix equation 

        Args:
            splitted_line (str): splitted line from netlist
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        node_a = int(splitted_line[1])
        node_b = int(splitted_line[2])
        Is = float(splitted_line[3])
        nVt = float(splitted_line[4])
        
        v0 = guess_v0diode
        G0 = Is*exp(v0/nVt)/nVt
        I0 = Is*(exp(v0/nVt) - 1) - G0*v0

        """print ("G0 =", G0)
        print ("I0 =", I0)
        input()"""

        Gn[node_a, node_a] +=  G0
        Gn[node_b, node_b] +=  G0
        Gn[node_a, node_b] += -G0
        Gn[node_b, node_a] += -G0

        In[node_a] += -I0
        In[node_b] +=  I0

        nodes = [node_a, node_b]
        return Gn, In, nodes

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
                or content[i][0] == "V" or content[i][0] == "D":

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

        self.n = len(nodes)
        return self.n, content

    def add_non_linear_components(self, var_content, Gn, In, v0):
        
        nodes_diodes = []
        i = 0
        for component in var_content:
            splitted_line = component.split(" ")
            
            if component[0] == "D":
                Gn, In, nodes = self.Add_Diode(splitted_line, Gn, In, v0[i])
                nodes_diodes.append(nodes)
            i += 1
                
        return Gn, In, nodes_diodes


    def add_components_time_variant(self, Gn, In, var_content, time):
        var_var_content = []
        for component in var_content:
            splitted_line = component.split(" ")
            
            if component[0] == "I" and splitted_line[3] == "SIN":
                In = self.Add_Isin(splitted_line, In, time)
            
            if component[0] == "V" and splitted_line[3] == "SIN":
                Gn, In = self.Add_Vsin(splitted_line, Gn, In, time)

            elif component[0] == "C":
                Gn, In = self.Add_C(splitted_line, Gn, In)

            elif component[0] == "K":
                Gn, In = self.Add_K(splitted_line, Gn, In)

            elif component[0] == "L":
                Gn, In = self.Add_L(splitted_line, Gn, In)
            
            elif component[0] == "D":
                var_var_content.append(component)
        
        return Gn, In, var_var_content

    def add_components_time_invariant(self, Gn, In, content):
        """Calls each "Add" function to add each component to the
        KCL equation. 

        Args:
            Gn (np.ndarray): Conductance matrix
            In (np.ndarray): Current vector
            content (list): Full content of netlist file in "readlines" format.

        Returns:
            np.ndarray, np.ndarray: new Gn and In matrix 
        """
        var_content = []
        for component in content:
            splitted_line = component.split(" ")
            
            if component[0] == "I":
                if splitted_line[3] == "DC":
                    In = self.Add_Idc(splitted_line, In)

            elif component[0] == "R":
                Gn = self.Add_R(splitted_line, Gn)
                
            elif component[0] == "G":
                Gn = self.Add_GmI(splitted_line, Gn)

            elif component[0] == "F":
                Gn, In = self.Add_F(splitted_line, Gn, In)
            
            elif component[0] == "E":
                Gn, In = self.Add_E(splitted_line, Gn, In)

            elif component[0] == "H":
                Gn, In = self.Add_H(splitted_line, Gn, In)
            
            elif component[0] == "V":
                if splitted_line[3] == "DC":
                    Gn, In = self.Add_VDC(splitted_line, Gn, In)
            
            # TIME VARIANTS #
            if component[0] == "I" and splitted_line[3] == "SIN":
                var_content.append(component)

            elif component[0] == "V" and splitted_line[3] == "SIN":
                var_content.append(component)

            elif component[0] == "C" or component[0] == "K" or\
            component[0] == "L" or component[0] == "D":
                var_content.append(component)

        return Gn, In, var_content

    def analyze(self, netlist, period, step, max_iter, tolerance):
        """Calculates the nodal voltages and variables of interest. 
        Stores values in attributes "e", "i_values", "vars_values".

        Args:
            netlist (str): Path for netlist file.

        Returns:
            np.ndarray: nodal voltages vector 
        """
        # Refresh all variables
        self.__init__()
        
        # Save delta(t) and period
        self.deltaT = step
        self.period = period

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

        # Add time invariant components
        invariant_Gn, invariant_In, var_content = self.add_components_time_invariant(Gn, In, content)
        
        for time in np.arange(0.0, period, step):
            self.vars = [] # Because add_components_time_variant adds the same vars each time

            # Add time variant components
            Gn_variant, In_variant, diode_components = self.add_components_time_variant(
                                                                invariant_Gn, 
                                                                invariant_In, 
                                                                var_content, 
                                                                time)

            ################################### ADDING DIODE ################################## 
            if time == 0.0:
                # Guess for diode "v0" 
                v0_diode = np.random.rand(len(diode_components))
                v0_diode = [0.5]
                # Set to 0.7 V if it is higher
                for i in range(len(v0_diode)):
                    if v0_diode[i] > 0.5:
                        v0_diode[i] = 0.5

            next_v0_diode = np.ones(len(diode_components)).flatten()
            epsilon = inf
            iters = 0
            while epsilon > tolerance and iters < max_iter: 
                Gn, In, nodes_diodes = self.add_non_linear_components(diode_components, Gn_variant, In_variant, v0_diode)

                # Removes node 0 (GND)
                Gn = Gn[1:, 1:]
                In = In[1:]
                
                # Solve the system Gn.e = In and get "e"
                vars_values = np.linalg.solve(Gn, In)
                vars_values = vars_values.flatten()

                epsilon = 0 # Assuming it to get the biggest epsilon between diodes

                for i in range(len(nodes_diodes)):
                    na = nodes_diodes[i][0]
                    nb = nodes_diodes[i][1]
                    
                    # Diode voltage from this iteration (Initial condition of the next)
                    next_v0_diode[i] = vars_values[na - 1] - vars_values[nb - 1]
                    
                    # Calculates difference between v0 and current v 
                    epsilon_test = abs(next_v0_diode[i] - v0_diode[i])

                    # Epsilon becomes the biggest epsilon among the diodes
                    if epsilon_test > epsilon:
                        epsilon = epsilon_test

                    # For next iteration v0 becomes the next v0     
                    v0_diode[i] = next_v0_diode[i]

                iters += 1
            
            # If it didn't converge
            if iters == max_iter:
                sys.exit("EXCEDEU NUMERO MAXIMO DE ITERAÇÔES !")
                
            ##############################################################################

            # Pick how many nodal voltage variables the system has  
            len_e = len(vars_values) - len(self.vars)

            # Creating list for all results of "e" through time. 
            if time == 0.0:
                for i in range(len_e):
                    self.e.append([vars_values[i]])
            else:
                for i in range(len_e):
                    self.e[i].append(vars_values[i])

            # Saves the "e" values in self.e0 (new initial condition)
            self.e0 = vars_values[:len_e] 
            
            # Saves i_values in self.i_values0 (new initial condition)
            self.i_values0 = vars_values[len_e:]

            # Saves all variables in self.vars_values0 (new initial condition)
            self.vars_values0 = vars_values[:]

            # Save i_values0 in dict format
            for i in range(len(self.i_values0)):
                self.i_vars0[self.vars[i]] = self.i_values0[i]

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
                print(f"e{i + 1} = {self.e[i]:.3f} V")
            
            for i in range(len(self.i_values)):
                print(f"{self.vars[i]} = {self.i_values[i]:.3f} A")

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

    def plot_nodes(self, plot_nodes): 
        x = np.arange(0.0, self.period, self.deltaT)
        fig, ax = plt.subplots(figsize=(5, 5))
        for node in plot_nodes:
            y = self.e[node - 1]
            
            color = (random.random(), random.random(), random.random())
            ax.plot(x, y, color=color, label=f'e{node}')
        
        ax.set_title("Nodal Voltages - Time domain", color='C0')
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Voltage (V)")
        ax.legend()
        plt.show()


def validate_list_of_plots_nodes(ctx, param, value):
    
    # Click library is executing this function twice, after returning 
    # some value. So, in the second time, "value" will be a list, and it
    # will return it right away.
    if isinstance(value, list):
        return value

    value_splitted = value.split(",")
    
    # Ignore "[]"
    try:
        if value_splitted[0][0] == "[":
            if len(value_splitted[0]) == 1:
                raise click.BadParameter("Invalid list. None value in first node", ctx, param)

            value_splitted[0] = value_splitted[0][1:]

        if value_splitted[-1][-1] == "]":
            if len(value_splitted[-1]) == 1:
                raise click.BadParameter("Invalid list. None value in last node", ctx, param)
            
            value_splitted[-1] = value_splitted[-1][:-1]
    except:
        raise click.BadParameter("Invalid list. It must have some value between \",\"", ctx, param)

    new_list_nodes = []
    for node in value_splitted:
        number = []
        for i in range(len(node)):      
            if node[i] == " ": continue
            if node[i] < "0" or node[i] > "9": # It means that it is not a number (ASCII)
                raise click.BadParameter("Invalid node. Nodes must be integers", ctx, param)

            number.append(node[i])
        if number == []:
            raise click.BadParameter("Invalid node. Node values must be integer and separeted by \",\"", ctx, param)

        new_list_nodes.append(int("".join(number)))

    if new_list_nodes == []:
        raise click.BadParameter("Invalid list. List is empty", ctx, param)

    # Verify if there are repeated nodes
    unique_nodes = set(new_list_nodes)
    if len(unique_nodes) != len(new_list_nodes):
        raise click.BadParameter("Invalid list. Can't have repeated nodes.", ctx, param)
    
    return new_list_nodes

@click.command()
@click.option("--netlist", help="Netlist file path", prompt=True, required=True)
@click.option("--period", help="Simulation time", prompt="Simulation time (s)", type=float)
@click.option("--step", help="Step time", prompt="Step time (s)", type=float)
@click.option("--max_iter", help="Maximum number of iterations in Newton-Raphson algorithm.", prompt=True, type=int)
@click.option("--epsilon_tolerance", help="Newton-Raphson tolerance for convergence", prompt=True, type=float)
@click.option("--plot_nodes", help="List of nodes to plot", prompt=True, type=click.UNPROCESSED, callback=validate_list_of_plots_nodes)
def main(netlist, period, step, max_iter, epsilon_tolerance, plot_nodes):
    """
    Calculates the nodal voltages of the circuit provided by the netlist.
    """
    
    period = period
    step = step

    ckt = Circuit()

    ckt.analyze(netlist, period, step, max_iter, epsilon_tolerance)
    
    ckt.plot_nodes(plot_nodes)
    
main()