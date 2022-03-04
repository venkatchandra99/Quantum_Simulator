# Quantum_Simulator
Classical Simulator for Quantum Computer using basic python.

## Creating object of Quantum Circuit class
qc = QuantumCircuit(n)

## Single qubit gates
qc.apply_gate('gate' , [qubit number])

## Two qubit gates
qc.apply_gate('two qubit gate' , [ [control bit] , [target qubit] ])

## Three qubit gates
qc.apply_gate('ccx' , [ [contro_1, control_2] , [target] ])

## Circuit to User defined gate
qc.to_gate(label = "name of gate")

## User defined gates
  qc2 = QuantumCircuit(m)   #where m >= n ,
  qc2.apply_gate("name of user defined gate" , [0 to n numbers])

## Getting Statevector
qc.statevector()

## Getting unitary matix
qc.get_unitary()

## Finding probabilities
qc.get_probabilities()

## Plotting bar graph
qc.draw_bargraph()
