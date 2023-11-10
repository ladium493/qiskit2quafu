# Introduction
This package converts the circuit of qiskit into a circuit supported by quafu
The converted circuit can be input to the cloud platform for processing

# Input and Output
## Input:
`qc: QuantumCircuit object of Qiskit`


`regName (optional): Modified register name`

`basis_gates (optional): The set of gates supported by the quafu chip`

`optimization_level (optional): optimization level`


## Output:
The function returns two objects, quafu_qc and qc_merge

`quafu_qc: The circuit instance of quafu`

`qc_merge: The converted qasm circuit in string form.`


# Example:
## Import packages
```python
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator
from qiskit2quafu import qiskit2quafu
from quafu import User, Task
```

## Create a circuit
```python
qiskitCircuit = QuantumCircuit(4)
qiskitCircuit.h([0,1,2,3])
qiskitCircuit.mcx([0,1,2],3)
qiskitCircuit.ry(0.6,0)
permute = Operator([[0, 0, 1, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 1, 0, 0, 0]])
qiskitCircuit.unitary(permute, [0,1,3], label='P')
qiskitCircuit.draw(output='mpl')
```

## Call the transfer function
```pythonfrom qiskit2quafu import qiskit2quafu
quafu_gates = ['cx','cy','cz','cp','u1','u2','u3','h','id','swap','cswap','p','rx','ry',
               'rz','x','y','z','s','sdg','t','tdg','sx','ccx','rxx','ryy','rzz']

quafuCircuit, transpiled_qasm = qiskit2quafu(qiskitCircuit)
quafuCircuit.draw_circuit()
```
## After transfer, the new circuit can be supported by quafu.
## To upload the task, refer to https://scq-cloud.github.io/


