using QuantumClifford, QuantumOpticsBase, CairoMakie

""" Questions to ask:
✓ 1. How to convert P"Z" to the large matrix view 
2. Why use mctrajectories for QEPOptimize circuit error rate check instead of petrajectories? or another trajectory?
3. brief explination of the meaning behind code distance of error correcting codes 


"""
### Representing pauli ops in efficient way
# Efficient simulation of pauli ops circuits
# 'condensed' representation
P"Z" # String macros
P"Z" * P"X"
P"ZZ" * P"XX"
P"Z" ⊗ P"X"
random_pauli(100) # Not an easy task !

# How is it efficient? Look inside:
bitstring.(UInt8.(P"XZ".xz)) # Get the .xz bit view of the pauli, remove some zeros at front, get bit view

bitstring.(UInt8.(P"XZIIYZ".xz)) 

# Since we use this a bunch, let's make a function
bits(arr) = bitstring.(UInt8.(arr))

bits(P"XZ".xz)

# How can we view commutativity?
comm(P"X",P"Z")
comm(P"X",P"X")

# So it is stored with UInt64... What about when it is too large?
a_billion = 1_000_000_000
l = random_pauli(a_billion)
# How?
l.xz 

r = random_pauli(a_billion)

# It's fast
using BenchmarkTools
@benchmark comm(l,r)
# Now he's showing off..
# 25 ms per comm op 
# seconds per qubits
25. / 1000 / 1e9 * 512

# how are these large ones stored?
random_pauli(67).xz

(P"X".xz,P"Z".xz) .|> bits
(P"X"*P"Z").xz |> bits
# and on the inside:
(P"X".xz .⊻ P"Z".xz) .|> bits
dump(P"X")
dump(P"iY")

### State representation
S"X" # <-- represents the one vector which is the eigenvector of the pauli X
# eigenstate which has the eigenvector +1
# it is a 'stabilizer' because it 'stabilizes' the pauli X
# using the Heisenberg representation 

# 'Ket' is from the QuantumOpticsBase
S"X" |> Ket # (|0⟩ + |1⟩)/√2
S"-X" |> Ket # (|0⟩ - |1⟩)/√2
S"XI II" |> Ket
S"XI IZ" 
S"ZI IZ" |> Ket # |00⟩
S"ZI -IZ" |> Ket # |10⟩
S"I" |> Ket # not really a valid system of equations.. This has many solutions 
S"XX XX" |> Ket # this is also not very useful, may solutions 
# We can still use these as 'mixed states'... will go over later.

# These are really just a matrix of binary numbers
S"XI IZ" |> stab_to_gf2

S"X" ⊗ S"Z" # Tensor product of two state vectors held by these constrains -> turns into block/diag tableaux

# When it is not seperable by block diag form, there is some sort of entanglement
S"XX ZZ" |> Ket # Bell state = |00⟩ + |11⟩
s = S"ZI IZ"
Ket(s)
s1 = S"ZZ IZ"
Ket(s1)
# Same state, why? they are system of linear equations. They are not uniquedly defined
# Same state, not the same tableaux
# We canonicalize with gaussian elimination
canonicalize!(s) == canonicalize!(s1)
s1 = S"ZZ IZ"
QuantumClifford.mul_left!(s1, 1, 2) # into row 1, multiply row 2
s1 = S"ZZ IZ"
s1[1]
s1[1] = s1[1]*s1[2]
s1[1]
# different ways to canonicalize

### Applying operations 
# Clifford ops map pauli -> pauli (preserve pauli ops)

tHadamard 
p = P"XYZ"
stab_to_gf2(p) # Decompose as the sum of binary basis vectors of pauli ops
stab_to_gf2(tab(tCNOT)) # extract tableaux representing cnot
tCNOT * S"XI IZ"
# And underneath:
gf2CNOT = stab_to_gf2(tab(tCNOT))'
# dump(tCNOT)
tCNOT

gf2State = stab_to_gf2(S"XI")'
gf2CNOT * gf2State
stab_to_gf2(tCNOT * S"XI")

# Converting to large representation
S"XI IZ" |> Ket
tHadamard |> Operator
tCNOT |> Operator 

# Applying unitaries
s = random_stabilizer(50)
# apply unitaries only on a few qubits 
# act only on 3rd qbut state
apply!(s, sHadamard(3))

# random_clifford(3) |> Operator


### Applying measurements 
project!(S"X",P"Z")[1]
project!(S"XI IZ",P"ZI")[1]
# we can randomize to act like a proper measurement 
apply!(S"XI IZ", sMZ(1))
# What about anticommuting with multiple rows?
project!(S"XI XZ",P"ZI")[1]

project!(S"X",P"X")[1] # State is the eigenstate of the measure op, so commutes and nothing happens (?)
project!(S"X",P"X")[3] #Sign 
project!(S"X",P"-X")[3] #Sign 

project!(S"XI IZ",P"XZ")[1] 


### Destabilizer formalism - improve performance, and help when intepreting the codes as error correcting
# We can 'precompile' some states and extra information 
S"XI IZ" |> typeof 

# Contains two sub-tableaux
S"XI IZ" |> MixedDestabilizer
# has 'destab' ops: each one anti-commutes with the corresponding stab 
# sort of like the encoding of anti-commuting ops in a basis 
# now, with this, commuting and non-commuting are equally fast (use this for simulation !!)
S"XI XI" |> canonicalize! # not really useful... But we can use it if we send it to a MixedDestabilizer
# GHZ state (|000⟩+|111⟩)
ghz_ = S"XXX ZZI IZZ"
S"XXX ZZI IZZ" |> Operator

apply!(MixedDestabilizer(ghz_), sMZ(1))

# Now le's try it again with a 'not pure' state
s = MixedDestabilizer(S"ZZI IZZ")

project!(s,P"XXX")[1]
# The projective measurement 'increased' the stabilizer --> initialzied into a more well-defined state 
using LinearAlgebra
project!(s,P"XXX")[1] |> rank # 3 independent constraints
MixedDestabilizer(S"ZZI IZZ") |> rank # only 2

traceout!(MixedDestabilizer(S"XXX ZZI IZZ"),1) # partial trace over the first qubit 
# lost a lot of rows in stab tableaux
# Because ghz is 'sensitive' destroying 1 qubit -> loose the ent. 
MixedDestabilizer(S"III") # Mixed state of all basis states


### Error correcting codes
using QuantumClifford.ECC 
# Check out QECCore
# This is a 9 qubit stabilizer 
p = parity_checks(Shor9()) # only 8 rows 
# These are just constraints on a Hilbert space 
# -> the constrained space is the code space 
p |> code_k # logical qubits
p |> code_n # physical qubits

naive_encoding_circuit(Shor9())
naive_syndrome_circuit(Shor9())

# We can build codes, but how can we bring them to practicality -> test their performance in real circuits? 
# estimating distance of code -> look at QECCore.distance


### Slightly non-clifford 
# can we do something like 1√2 S"X" + 1/√2im S"Z" ?
s = GeneralizedStabilizer(S"X")
pcT # pauli channel, T gate 
pcT |> Operator
normalize!(pcT |> Operator)
apply!(s,pcT) 
apply!(s,pcT) |> Operator

