from sys import argv

assert(len(argv) > 4)
q_coeffs_b = dict(eval(argv[1]))
codeword_b = list(eval(argv[2]))
beta_coeff_list = list(eval(argv[3]))
irred_coeff_list = list(eval(argv[4]))


irred = 0
for i,coeff in enumerate(irred_coeff_list):
	irred += coeff * x^i
k.<x> = GF(2^(len(irred_coeff_list)-1),'x', irred)
R.<y,z> = k[]
beta = 0
for i,coeff in enumerate(beta_coeff_list):
	beta += coeff * x^i
assert(beta in k)

Q = 0
degy = 0
for i,j in q_coeffs_b.keys():
	degy = max(degy,i)
	if q_coeffs_b[(i,j)] != 0:
		Q += (beta^q_coeffs_b[(i,j)]) * y^i * z^j

assert(Q in R)

for i in range(len(codeword_b)):
	if codeword_b[i] == 0:
		assert(Q(beta^i,0) == 0)
	else:
		assert(Q(beta^i,beta^codeword_b[i]) == 0)

print("Got Polynomial ",Q, file=sys.stderr)
Q_factors = Q.factor()



msg_list = []
print("Q has the following factors: ",file=sys.stderr)
for fac,expo in Q_factors:
	print(fac,file=sys.stderr)
	if z in fac.variables() and z not in (fac-z).variables():
		msg_list.append(-(fac-z))

msg_list_b = []
for m in msg_list:
	m_b = [0 for i in range(degy+1)]
	for i,mono in enumerate(m.monomials()):
		for p in range(degy+1):
			if(y^p == mono):
				tmp = m.coefficients()[i]
				
				for p2 in range(1,1<<(len(irred_coeff_list)-1)):
					if beta^p2 == tmp:
						m_b[p] = p2

	msg_list_b.append(m_b)


print(msg_list_b)