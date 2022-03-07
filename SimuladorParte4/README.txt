Instruções de uso

Bibliotecas necessárias:
    - numpy
    - click
    - matplotlib

Como executar:
    python simulador_parte4.py --netlist=<path da netlist> --period=<float> 
    --step=<float> --max_iter=<int> --epsilon_tolerance=<float> --plot_nodes=<list>

Exemplo:
    python simpulador_parte3.py --netlist=netlist0.txt --period=3 --step=5e-4 
    --max_iter=1000 --epsilon_tolerance=1e-4 --plot_nodes=1,2,3

Ou executar o programa sem passar argumentos.

Exemplo:
    python simulador_parte4.py
    
    >> Netlist: netlist0.txt
    >> Simulation time (s): 3
    >> Step time (s): 5e-4
    >> Max iter: 5e-4
    >> Epsilon tolerance: 5e-4
    >> Plot nodes: 1,2,3     OU    >> Plot nodes: [1,2,3]