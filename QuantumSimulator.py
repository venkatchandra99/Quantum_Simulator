import matplotlib.pyplot as plt
import math

class QuantumCircuit():
    def __init__(self,num_of_Qubits,num_of_Cbits=0):
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
        self.S=[[1,0],[0,complex(i)]]
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
    
    def apply_gate(self,gate,index,theta=0):
        if gate=='x' or gate=='cx' or gate =='ccx':
            self.gates_executed.append([index,gate])
            gate=self.X
        elif gate=='h'or gate=='ch':
            self.gates_executed.append([index,gate])
            gate=self.H
        elif gate=='s' or gate=='cs':
            self.gates_executed.append([index,gate])
            gate=self.S
        elif gate=='t' or gate=='ct':
            self.gates_executed.append([index,gate])
            gate=self.t
        elif gate=='y' or gate=='cy':
            self.gates_executed.append([index,gate])
            gate=self.Y
        elif gate=='z' or gate=='cz':
            self.gates_executed.append([index,gate])
            gate=self.Z
        elif gate=='rx' or gate=='crx':
            self.gates_executed.append([index,gate,"theta: "+str(theta*180//math.pi)])
            gate=[[math.cos(theta/2),complex(0,-math.sin(theta/2))],[complex(0,-math.sin(theta/2)),math.cos(theta/2)]]
        elif gate=='rz' or gate=='crz':
            self.gates_executed.append([index,gate,"theta: "+str(theta*180//math.pi)])
            gate=[[complex(math.cos(theta/2),-math.sin(theta/2)),0],[0,complex(math.cos(theta/2),math.sin(theta/2))]]
        elif gate=='ry' or gate=='cry':
            self.gates_executed.append([index,gate,"theta: "+str(theta*180//math.pi)])
            gate=[[math.cos(theta/2),-math.sin(theta/2)],[math.sin(theta/2),math.cos(theta/2)]]
        elif gate=='swap':
            self.apply_gate('cx',index,theta=0)
            self.apply_gate('cx',index[::-1],theta=0)
            self.apply_gate('cx',index,theta=0)
            return 0
            
        else:
            raise Exception("Specify valid gate")
        
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

        # Two Qubit gates 
        elif len(index)==2:
            control=index[0]
            target=index[1][0]
            for ctrl in control:
                if (ctrl > self.num_of_Qubits - 1) or (target > self.num_of_Qubits - 1):
                    raise Exception("Qubit index is greater than Number of Qubits")
    
            n_qubits=self.num_of_Qubits
            res=[]
            basis_states=[]
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
    
    def get_unitary(self):
        return self.unitary
    
    def get_probabilities(self):
        prob={}
        _state_vector=self.matrix_mul(self.unitary,self.state_vector)
        for i in range(len(self.basis)):
            prob[self.basis[i]]=abs(_state_vector[i][0]**2)
        return prob
            
    def draw_bargraph(self):
        prob = self.get_probabilities()
        x=list(prob.keys())
        y=list(prob.values())
        plt.xlabel('Basis State')
        plt.ylabel('Probability')
        if len(self.gates_executed) < 9:
            plt.title("Gates: "+str(self.gates_executed))
        else:
            plt.title("Circuit Probability Graph")
        plt.bar(x,y)
        plt.show()
    



qc=QuantumCircuit(3)
qc.apply_gate('x',[0])
qc.apply_gate('x',[1])
qc.apply_gate('rx',[1],math.pi/2)
qc.apply_gate('cx',[[0],[2]])
qc.apply_gate('ccx',[[0,1],[2]])
qc.apply_gate('swap',[[1],[0]])

print(qc.get_unitary())
print(qc.statevector())
print(qc.get_probabilities())
qc.draw_bargraph()
