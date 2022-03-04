import matplotlib.pyplot as plt
import math

global userdefined_gates
userdefined_gates={}

class Gates():
    def __init__(self,theta=0):
        self._gates={}

        #basic gates
        self._gates['I']=[[1,0],[0,1]]
        self._gates[('x','cx','ccx')]=[[0,1],[1,0]]
        self._gates[('y','cy')]=[[0,complex(0,-1)],[complex(0,1),0]]
        self._gates[('z','cz')]=[[1,0],[0,-1]]
        self._gates[('s','cs')]=[[1,0],[0,complex(0,1)]]
        self._gates[('h','ch')]=[[1/math.sqrt(2),1/math.sqrt(2)],[1/math.sqrt(2),-1/math.sqrt(2)]]
        self._gates[('t','ct')]=[[1,0],[0,math.cos(math.pi/4)+complex(0,math.sin(math.pi/4))]]

    #returns gates
    def gate(self,g,theta=0):
        if g in ['rx','crx']:
            return [[math.cos(theta/2),complex(0,-math.sin(theta/2))],[complex(0,-math.sin(theta/2)),math.cos(theta/2)]]
        elif g == ['ry','cry']:
            return [[math.cos(theta/2),-math.sin(theta/2)],[math.sin(theta/2),math.cos(theta/2)]]
        elif g == ['rz','crz']:
            return [[complex(math.cos(theta/2),-math.sin(theta/2)),0],[0,complex(math.cos(theta/2),math.sin(theta/2))]]
        elif g in self._gates.keys():
            return self._gates[g]
        else:
            for gates in self._gates.keys():
                if isinstance(gates,tuple):
                    if g in gates:
                        return self._gates[gates]
            if g in userdefined_gates.keys():
                return userdefined_gates[g]
                    
    def get_basic_gates(self):
        return ['I','x','cx','ccx','y','cy','z','cz','s','cs','h','ch','t','ct','rz','ry','rx','crx','crz','cry']
    def get_userdefined_gates(self):
        return userdefined_gates.keys()
    def get_allgates(self):
        all_gates = self.get_basic_gates()+list(userdefined_gates.keys())
        return all_gates
        
    #adds userdefined gates                
    def append_gate(self,label,unitary_matrix):
        userdefined_gates[label]=unitary_matrix
        
                

class QuantumCircuit(Gates):
    def __init__(self,num_of_Qubits,num_of_Cbits=0):
        Gates.__init__(self)
        self.num_of_Qubits=num_of_Qubits
        self.num_of_Cbits=num_of_Cbits
        self.gates_executed=[]
        #basis_states
        self.basis=[]
        
        for i in range(2**self.num_of_Qubits):
            self.basis.append(self.decimal_to_binary(i,self.num_of_Qubits))
        
        #initial state -> |0>
        self.zero_state=[[1],[0]]

        #basic gates
        self.I=[[1,0],[0,1]]
        self.X=[[0,1],[1,0]]
        self.Y=[[0,complex(0,-1)],[complex(0,1),0]]
        self.Z=[[1,0],[0,-1]]
        self.S=[[1,0],[0,complex(0,1)]]
        self.H=[[1/math.sqrt(2),1/math.sqrt(2)],[1/math.sqrt(2),-1/math.sqrt(2)]]
        self.T=[[1,0],[0,math.cos(math.pi/4)+complex(0,math.sin(math.pi/4))]]

        self.unitary=[[1]]
        self.state_vector=[[1]]
        for i in range(self.num_of_Qubits):
            self.unitary=self.tensor_product(self.unitary,self.I)
            self.state_vector=self.tensor_product(self.state_vector,self.zero_state)
        
    #Converts decimal to binary 
    def decimal_to_binary(self,num,i):
        return format(num, f"0{i}b")

    #Converts binary to decimal
    def binary_to_decimal(self,num):
        if isinstance(num,int):
            num=str(num)
        return int(num,2)

    #Calculates tensor product of given matrices or vectors(by default considered as coloumn vector)
    def tensor_product(self,A,B,row_A=False,row_B=False):
        if isinstance(A[0],int):
            if row_A:
                A=[A]
            else:
                for i in range(len(A)):
                    A[i]=[A[i]]
        if isinstance(B[0],int):
            if row_B:
                B=[B]
            else:
                for i in range(len(B)):
                    B[i]=[B[i]]
        
        res=[]
        for i in range(len(A)):
            for k in range(len(B)):
                res.append([])
                for j in range(len(A[i])):
                    for l in range(len(B[k])):
                        res[len(res)-1].append(A[i][j]*B[k][l])
        return res

    #Calculates product of two matrices
    def matrix_mul(self,mat1,mat2):
        res = [[0 for x in range(len(mat2[0]))] for y in range(len(mat1))]  
        for i in range(len(mat1)): 
            for j in range(len(mat2[0])): 
                for k in range(len(mat2)): 
                    res[i][j] += mat1[i][k] * mat2[k][j] 
        return res

    def apply_gate(self,_gate,index,theta=0):
        if _gate=='swap':
            self.apply_gate('cx',index,theta=0)
            self.apply_gate('cx',index[::-1],theta=0)
            self.apply_gate('cx',index,theta=0)
            return 
        else:
            gate = super().gate(_gate,theta)

        if gate is None:
            raise Exception("Specify valid gate")
        else:
            if _gate in ['rz','ry','rx','crx','crz','cry']:
                self.gates_executed.append([index,_gate+" "+str(round(theta,3))])
            else:
                self.gates_executed.append([index,_gate])

        #for basic gates
        if _gate in super().get_basic_gates():
            
            #Single Qubit gates
            if len(index)==1:
                q_index = index[0]
                if q_index > self.num_of_Qubits - 1:
                    raise Exception("Qubit index is greater than Number of Qubits")
            
                t=[[1]]
                for i in range(self.num_of_Qubits-1,-1,-1):
                    if i==q_index:
                        t=self.tensor_product(t,gate)
                    else:
                        t=self.tensor_product(t,self.I)
                self.unitary = self.matrix_mul(t,self.unitary)

            # Multi Qubit gates 
            elif len(index)==2:
                control=index[0]
                target=index[1][0]
                for ctrl in control:
                    if (ctrl > self.num_of_Qubits - 1) or (target > self.num_of_Qubits - 1):
                        raise Exception("Qubit index is greater than Number of Qubits")
        
                n_qubits=self.num_of_Qubits
                res=[]
                basis_states=[]

                #Applying the the specified gate on target qubit with respect to the control qubit
                for i in range(2**n_qubits):
                    res.append([])
                    deci_before = self.decimal_to_binary(i,n_qubits)
                    flag = 1
                    for ctrl in control:
                        if deci_before[len(deci_before)-1-ctrl]!='1':
                            flag=0
                    if flag == 1:
                        if deci_before[len(deci_before)-1-target]=='0':
                            deci_after=deci_before[0:len(deci_before)-1-target]+"1"+deci_before[len(deci_before)-target:]
                            int_after = self.binary_to_decimal(deci_after)
                            for k in range(2**n_qubits):
                                if k == i:
                                    res[i].append(gate[0][0])
                                elif k == int_after:
                                    res[i].append(gate[0][1])
                                else:
                                    res[i].append(0)
                        else:
                            deci_after=deci_before[0:len(deci_before)-1-target]+"0"+deci_before[len(deci_before)-target:]
                            int_after = self.binary_to_decimal(deci_after)
                            for k in range(2**n_qubits):
                                if k == int_after:
                                    res[i].append(gate[1][0])
                                elif k == i:
                                    res[i].append(gate[1][1])   
                                else:
                                    res[i].append(0)            
                    else:
                        for k in range(2**n_qubits):
                            if k == i:
                                res[i].append(1)
                            else:
                                res[i].append(0)
                self.unitary = self.matrix_mul(res,self.unitary)

        #for user defined gates
        else:
            for num in range((len(gate)//2)+1):
                if (2**num) == len(gate):
                    if self.num_of_Qubits < num:
                        raise Exception(f"{_gate} gate requires atleat {num} Qubits")
                    else:
                        break
            if len(index) != num:
                raise Exception(f"{_gate} gate works on {num} Qbits")
            
            #Creating an Identity matrix of size N  (N=2**num_of_qubits)
            matrix, gate_basis = [], []
            for i in range(2**self.num_of_Qubits):
                matrix.append([])
                for j in range(2**self.num_of_Qubits):
                    if i==j:
                        matrix[i].append(1)
                    else:
                        matrix[i].append(0)

            #getting bais states in according to the no.of qubits the userdefined gate works on.
            for i in range(2**len(index)):
                gate_basis.append(self.decimal_to_binary(i,len(index)))

            #Fitting the MxM userdefined gate in NxN identity matrix (N<=M)
            for row_i in range(len(self.basis)):
                tem_row = ""
                tem_row_2=""
                for k in range(len(self.basis[row_i])):
                    if k in index:
                        tem_row=tem_row + self.basis[row_i][::-1][k]
                    else:
                        tem_row_2=tem_row_2 + self.basis[row_i][::-1][k]
                tem_row=tem_row[::-1]
                for row_j in gate_basis:
                    if row_j == tem_row: 
                        matrix_row = self.binary_to_decimal(self.basis[row_i])
                        gate_row = self.binary_to_decimal(row_j)
                        for col_i in range(len(self.basis)):
                            tem_col = ""
                            tem_col_2 = ""
                            for k in range(len(self.basis[col_i])):
                                if k in index:
                                    tem_col=tem_col + self.basis[col_i][::-1][k]
                                else:
                                    tem_col_2=tem_col_2 + self.basis[col_i][::-1][k]
                                    
                            tem_col=tem_col[::-1]
                            for col_j in gate_basis:
                                if (col_j == tem_col) and (tem_row_2 == tem_col_2):
                                    matrix_col = self.binary_to_decimal(self.basis[col_i])
                                    gate_col = self.binary_to_decimal(col_j)
                                    matrix[matrix_row][matrix_col] = gate[gate_row][gate_col]
                                    break
                        break
            self.unitary = self.matrix_mul(matrix,self.unitary)

    #returns Statevector
    def statevector(self):
        _state_vector=self.matrix_mul(self.unitary,self.state_vector)
        for col in range(len(_state_vector)):
            for prob in range(len(_state_vector[col])):
                try:
                    _state_vector[col][prob]=round(_state_vector[col][prob],15)
                except TypeError as e:
                    _state_vector[col][prob]=round(_state_vector[col][prob].real, 15) + round(self.state_vector[col][prob].imag, 15) * 1j
##                    print(e)
        return _state_vector

    #returns unitary matrix
    def get_unitary(self):
        return self.unitary

    #returns probabilities
    def get_probabilities(self):
        prob={}
        _state_vector=self.matrix_mul(self.unitary,self.state_vector)
        for i in range(len(self.basis)):
            prob[self.basis[i]]=abs(_state_vector[i][0]**2)
        return prob

    #plots graph of probabilities between basis states
    def draw_bargraph(self):
        prob = self.get_probabilities()
        x=list(prob.keys())
        y=list(prob.values())
        temp=len(y)
        while temp:
            if y[temp-1]==0:
                del y[temp-1]
                del x[temp-1]
            temp=temp-1
        plt.xlabel('Basis State')
        plt.ylabel('Probability')
        if len(self.gates_executed) < 9:
            plt.title("Gates: "+str(self.gates_executed))
        else:
            plt.title("Circuit Probability Graph")
        plt.bar(x,y)
        plt.show()

    #converts circuit to gate
    #circuit's unitary matrix is added to gates
    def to_gate(self,label):
        if (label in super().get_allgates()) or label == 'swap':
            print("Gate already exists")
        else:
            unitary_matrix=self.get_unitary()
            super().append_gate(label,unitary_matrix)
        

#Example
#Implementing Cuccaro adder circuit

#Majority Circuit
maj_c = QuantumCircuit(3)
maj_c.apply_gate('cx',[[2],[1]])
maj_c.apply_gate('cx',[[2],[0]])
maj_c.apply_gate('ccx',[[0,1],[2]])
maj_c.to_gate(label="MAJ")

#UnMajority Circuit
uma3_c = QuantumCircuit(3)
uma3_c.apply_gate('x',[1])
uma3_c.apply_gate('cx',[[0],[1]])
uma3_c.apply_gate('ccx',[[0,1],[2]])
uma3_c.apply_gate('x',[1])
uma3_c.apply_gate('cx',[[2],[0]])
uma3_c.apply_gate('cx',[[2],[1]])
##uma3_c.apply_gate('cx',[[0],[1]])
uma3_c.to_gate(label="UMA3")

#2-Qubit Cuccar adder circuit
C_adder = QuantumCircuit(6)
C_adder.apply_gate('MAJ',[0,1,2])
C_adder.apply_gate('MAJ',[2,3,4])
C_adder.apply_gate('cx',[[4],[5]])
C_adder.apply_gate('UMA3',[2,3,4])
C_adder.apply_gate('UMA3',[0,1,2])
C_adder.to_gate(label="C_adder(2)")

#adding 3(11) and 1(01)
Add = QuantumCircuit(6)

#3(11)
Add.apply_gate('x',[2])
Add.apply_gate('x',[4])

#1(01)
Add.apply_gate('I',[3])
Add.apply_gate('x',[1])


Add.apply_gate('C_adder(2)',[0,1,2,3,4,5])
prob = Add.get_probabilities()

x,y =[],[]
for i in prob.keys():
    if prob[i]!=0:
        j=i[::-1]
        val = j[1]+j[3]+j[5]
        val = val[::-1]
        x.append(val)
        y.append(prob[i])

plt.xlabel('Basis State ')
plt.ylabel('Probability')
print(x)
plt.title("2Qubit cucarro adder --> 11 (3) + 01 (1) = 100 (4)")
plt.bar(x,y)
plt.show()
##Add.draw_bargraph()


