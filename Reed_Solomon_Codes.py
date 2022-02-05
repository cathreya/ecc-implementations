from sympy import *
import numpy as np
from math import log
from pandas import *
from math import ceil
x, y, z = symbols('x,y,z')
# init_printing(use_unicode=False, wrap_line=False)

print(
"""
#################################################################
# Q1 (a) Finding an irreducible polynomial of degree 8 over F_2 #
#################################################################
"""
)

def get_irreducible(deg = 8):
	"""
	get_irreducible brute forces over all polynomials of 
	degree deg and attempts to factorize.	
	"""
	for i in range(1,1<<(deg-1)):
		f = 0
		for j in range(deg-1):
			if((i&(1<<j))):
				f += x**j

		f += x**deg
		f = f.as_poly(domain=FF(2))
		if f.is_irreducible:
			return f


def is_primitive(f, irred):
	"""
	is_primitive returns true if f generates all 
	non-zero elements in F_2[x]/irred
	"""
	freq = {}
	cur = f
	cnt = 0
	deg = irred.degree()
	for i in range(1<<deg):
		coeffs = tuple(cur.all_coeffs())
		if coeffs not in freq.keys():
			cnt += 1
		freq[coeffs] = 1
		
		cur *= f
		_,cur = div(cur,irred, domain=FF(2))

	return cnt == (1<<deg)-1


def find_primitive_element(irred):
	"""
	find_primitive_element brute forces over polynomials
	and checks whether they're primitive for F_2[x]/irred
	"""
	deg = irred.degree()
	for i in range(1,1<<deg):
		f = Poly(0,x, domain = FF(2))
		for j in range(deg):
			if((i&(1<<j))):
				f += x**j

		f = f.as_poly(domain=FF(2))

		if is_primitive(f, irred):
			return f


def fast_pow(cur, k, irred):
	"""
	computes cur^k mod irred in O(log k)
	"""
	if k == 0:
		return Poly(1,x,domain=FF(2))
	half = fast_pow(cur,k//2,irred)
	if k%2 == 0:
		res = half*half
		_,res = div(res,irred,domain=FF(2))
		return res
	else:
		res =  cur*half*half
		_,res = div(res,irred,domain=FF(2))
		return res
		
f = get_irreducible()

print(f," is an irreducible polynomial of degree 8 in F_2[x]")

print(
"""
#################################################################
# Q1 (b) Finding a primitive element                            #
#################################################################
"""
)

beta = find_primitive_element(f)
beta34 = fast_pow(beta, 34, f)
beta20 = fast_pow(beta,20,f)
beta54prod = beta20*beta34
_,beta54prod = div(beta54prod, f, domain=FF(2))

beta54 = fast_pow(beta,54,f)
_,beta54 = div(beta54, f, domain=FF(2))

assert(beta54 == beta54prod)
print(beta54," == ", beta54prod)

print(
"""
#################################################################
# Q2 Calculating encoded RS code vectors                        #
#################################################################
"""
)


def gen_vandermonde(n,k, beta, irred):
	"""
	Outputs an nxk vandermonde matrix
	"""
	M = []
	for i in range(n):
		row = []
		cur = fast_pow(beta,i,irred)
		for j in range(k):
			row.append(fast_pow(cur,j,irred))
		M.append(row)

	return M

def mat_vec_prod(M,v,irred):
	"""
	Multiply a matrix and vector of polynomials 
	"""
	assert(len(M[0]) == len(v))
	res = [0 for i in range(len(M))]

	for i in range(len(M)):
		for j in range(len(v)):
			res[i] += M[i][j]*v[j]
			_,res[i] = div(res[i],irred,domain=FF(2))
	return res

print(
"""
###################
# (a) [n=10, k=4] #
###################
"""
)

# For n = 10 we require field size 2^m > 10, let m = 4
f = get_irreducible(4)
beta = find_primitive_element(f)

G = gen_vandermonde(10,4,beta,f)
message = [
	Poly(1,x,domain=FF(2)),
	fast_pow(beta,2,f),
	Poly(0,x,domain=FF(2)),
	Poly(0,x,domain=FF(2)),
]
codeword = mat_vec_prod(G,message,f)

print(np.array(codeword)[:,None])

print(
"""
####################
# (b) [n=20, k=10] #
####################
"""
)
# For n = 10 we require field size 2^m > 20, let m = 5
f = get_irreducible(5)
beta = find_primitive_element(f)

G = gen_vandermonde(20,10,beta,f)
message = [
	fast_pow(beta,4,f),
	Poly(0,x,domain=FF(2)),
	Poly(0,x,domain=FF(2)),
	fast_pow(beta,5,f),
	Poly(0,x,domain=FF(2)),
	Poly(0,x,domain=FF(2)),
	Poly(0,x,domain=FF(2)),
	Poly(0,x,domain=FF(2)),
	Poly(0,x,domain=FF(2)),
	fast_pow(beta,10,f),
]
codeword = mat_vec_prod(G,message,f)
print(np.array(codeword)[:,None])

print(
"""
####################
# (c) [n=50, k=25] #
####################
"""
)
# For n = 50 we require field size 2^m > 50, let m = 6
f = get_irreducible(6)
beta = find_primitive_element(f)


G = gen_vandermonde(50,25,beta,f)
message = [	Poly(0,x,domain=FF(2)) for i in range(25)]
message[0] = fast_pow(beta,20,f)
message[8] = fast_pow(beta,4,f)
message[12] = fast_pow(beta,44,f)
message[14] = fast_pow(beta,30,f)
message[18] = fast_pow(beta,15,f)
message[23] = fast_pow(beta,3,f)

codeword = mat_vec_prod(G,message,f)
print(np.array(codeword)[:,None])


print(
"""
#################################################################
# Q3 RS decoding using the Welch Berlekamp Algorithm            #
#################################################################
"""
)

n = int(input("Enter n "))
k = int(input("Enter k "))
siz = int(log(n,2))+1
f = get_irreducible(siz)
beta = find_primitive_element(f)

message = [	Poly(0,x,domain=FF(2)) for i in range(k)]

non_zero = int(input("Enter the number of non-zero coefficients "))
for i in range(non_zero):
	j = int(input("Enter power of X "))
	assert(j<len(message))
	v = int(input("Enter p, where coefficient in the form beta^p (beta is the primitive element) "))
	message[j] = fast_pow(beta,v,f)

print(
"""
####################
# Message Vector   #
####################
"""
)
print(np.array(message)[:,None])

G = gen_vandermonde(n,k,beta,f)
codeword = mat_vec_prod(G,message,f)
print(
"""
####################
# Codeword         #
####################
"""
)
print(np.array(codeword)[:,None])

dmin = n-k+1
t = (dmin-1)//2
e = int(input("Enter the number of errors(at most {}): ".format(t)))
assert(e <= t)
for i in range(e):
	j = int(input("Enter position "))
	assert(j<len(codeword))
	v = int(input("Enter p, used to derive a polynomial from F_2[x] of the form beta^p (beta is the primitive element) "))
	codeword[j] = fast_pow(beta,v,f)

print(
"""
####################
# Received Vector  #
####################
"""
)

print(np.array(codeword)[:,None])


# deg(E) = t
# deg(N) = t+k-1
G1 = gen_vandermonde(n,t+1,beta,f)
G2 = gen_vandermonde(n,t+k,beta,f)

constraints = []
for i in range(n):
	row = []
	for j in range(t):
		res = G1[i][j]*codeword[i]
		_,res = div(res,f,domain=FF(2))
		row.append(res)

	for j in range(t+k):
		res = f-G2[i][j]
		_,res = div(res,f,domain=FF(2))
		row.append(res)

	res = G1[i][t]*codeword[i]
	_,res = div(res,f,domain=FF(2))
	res = f-res
	_,res = div(res,f,domain=FF(2))
	row.append(res)
	constraints.append(row)


def inv(a, irr_poly):
	# print("trying to invert", a)
	return invert(a,irr_poly)

def print_mat(a):
	copy = []
	for i in range(len(a)):
		copy_row=[]
		for j in range(len(a[0])):
			copy_row.append(a[i][j].as_expr())

		copy.append(copy_row)
	print(DataFrame(copy))

def GaussElim(A, irred, silent = False):
	"""
	Solves the system of linear equations given by A 
	(constant specified by the last column)
	over polynomials in F_2[x] mod irred
	"""

	mat = A
	row = [-1 for i in range(len(A[0]))]
	R = len(mat)
	C = len(mat[0])-1
	r = 0
	for c in range(C):
		if not silent:
			print("Eliminating Column {}\n".format(c))
			print_mat(mat)
			print("\n")
		k = r
		while k<R and mat[k][c] == Poly(0,x,domain=FF(2)):
			k+=1
		if k==R:
			if not silent:
				print("No pivot found")
			continue

		if not silent:
			print("Pivot: ", k, mat[k][c].as_expr())
		mat[k],mat[r] = mat[r],mat[k]

		mod_inv = inv(mat[r][c],irred)
		for i in range(R):
			if i != r:
				w = mat[i][c]*mod_inv
				_,w = div(w, irred, domain=FF(2))
				w = irred - w
				_,w = div(w, irred, domain=FF(2))

				for j in range(C+1):
					mat[i][j] = (mat[i][j] + mat[r][j] * w)
					_,mat[i][j] = div(mat[i][j],irred,domain=FF(2))

		row[c] = r
		r += 1


	if not silent:
		print("Final Constraint Matrix\n")
		print_mat(mat),
		print("\n")
	ans = [0 for i in range(C)]
	for i in range(C):
		r = row[i]
		if r == -1:
			ans[i] = Poly(0,x,domain=FF(2))
		else:
			ans[i] = (mat[r][C] * inv(mat[r][i],irred))
			_,ans[i] = div(ans[i],irred,domain=FF(2))

	return ans,row


eandn,_ = GaussElim(constraints,f)

ex = [eandn[i] for i in range(t)]
ex.append(Poly(1,x,domain=FF(2)))
nx = [eandn[i] for i in range(t,len(eandn))]

print(
"""
####################
# Error Polynomial #
####################
"""
)
print(np.array(ex)[:,None])

print(
"""
####################
# n(x) Polynomial  #
####################
"""
)
print(np.array(nx)[:,None])


def poly_div(u,v,irred):
	"""
	Returns q,r s.t. v = qu + r
	where u and v are polynomials in F_2[x] mod irred
	"""
	i = len(u)-1
	j = len(v)-1
	mx = [Poly(0,x,domain=FF(2)) for i in range(len(u)-len(v)+1)]
	while(j>=0 and v[j] == Poly(0,x,domain=FF(2))):
		j-=1
	if(j == 0):
		return 0,0
	while(i>=j):
		if(u[i] == Poly(0,x,domain=FF(2))):
			i-=1
			continue

		tmp = u[i]*inv(v[j],irred)
		_,tmp = div(tmp,irred,domain=FF(2))
		mx[i-j] = tmp

		for p in range(len(v)):
			tmp2 = tmp*v[p]
			_,tmp2 = div(tmp2,irred,domain=FF(2))
			tmp2 = irred - tmp2
			_,tmp2 = div(tmp2,irred,domain=FF(2))
			u[p+i-j] += tmp2
			_,u[p+i-j] = div(u[p+i-j],irred,domain=FF(2))

		i-=1

	return mx,u 

mx,rem = poly_div(nx,ex,f)
if mx ==0 and rem == 0:
	pass
	# need to handle the case where E = 0


print(
"""
####################
# Decoded Message  #
####################
"""
)
print(np.array(mx)[:,None])

print(
"""
####################
# Original Message #
####################
"""
)
print(np.array(message)[:,None])



print(
"""
#################################################################
# Q4 List Decoding Using Algorithm 1                            #
#################################################################
"""
)

n = int(input("Enter n "))
k = int(input("Enter k "))
siz = int(log(n,2))+1
f = get_irreducible(siz)
beta = find_primitive_element(f)

message = [	Poly(0,x,domain=FF(2)) for i in range(k)]

non_zero = int(input("Enter the number of non-zero coefficients "))
for i in range(non_zero):
	j = int(input("Enter power of X "))
	assert(j<len(message))
	v = int(input("Enter p, where coefficient in the form beta^p (beta is the primitive element) "))
	message[j] = fast_pow(beta,v,f)

print(
"""
####################
# Message Vector   #
####################
"""
)
print(np.array(message)[:,None])

G = gen_vandermonde(n,k,beta,f)
codeword = mat_vec_prod(G,message,f)
print(
"""
####################
# Codeword         #
####################
"""
)
print(np.array(codeword)[:,None])

dmin = n-k+1
t = (dmin-1)//2
e = int(input("Enter the number of errors (upto {}+1): ".format(t)))
assert(e <= t+1)
for i in range(e):
	j = int(input("Enter position "))
	assert(j<len(codeword))
	v = int(input("Enter p, used to derive a polynomial from F_2[x] of the form beta^p (beta is the primitive element) "))
	codeword[j] = fast_pow(beta,v,f)

print(
"""
####################
# Received Vector  #
####################
"""
)

print(np.array(codeword)[:,None])



degx = ceil(sqrt(n*(k-1)))
# degy = ceil(n/degx)
degy = 2

constraints = []

for row in range(n):
	alpha = fast_pow(beta,row,f)
	y = codeword[row]
	c_row = []
	for i in range(degx+1):
		for j in range(degy + 1):
			tmp = fast_pow(alpha,i,f) * fast_pow(y,j,f)
			_,tmp = div(tmp,f,domain=FF(2))
			c_row.append(tmp)
	c_row.append(Poly(0,x,domain=FF(2)))
	constraints.append(c_row)


q_coeffs,free = GaussElim(constraints, f, silent=True)

freei = -1
freej = -1
for i in range(degx,-1,-1):
	fg = 0
	for j in range(degy,-1,-1):
		if free[i*(degy+1)+j] == -1:
			freei = i
			freej = j
			fg = 1
			break
	if fg:
		break

constraints = []
for row in range(n):
	alpha = fast_pow(beta,row,f)
	y = codeword[row]
	c_row = []
	for i in range(degx+1):
		for j in range(degy + 1):
			if i == freei and j == freej: 
				continue
			tmp = fast_pow(alpha,i,f) * fast_pow(y,j,f)
			_,tmp = div(tmp,f,domain=FF(2))
			c_row.append(tmp)

	tmp = fast_pow(alpha,freei,f) * fast_pow(y,freej,f)
	_,tmp = div(tmp,f,domain=FF(2))
	tmp = f-tmp
	_,tmp = div(tmp,f,domain=FF(2))
	c_row.append(tmp)
	constraints.append(c_row)

q_coeffs,_ = GaussElim(constraints, f)
q_coeffs.insert(freei*(degy+1)+freej, Poly(1,x,domain=FF(2)))


q_coeffs_b = {}

for i in range(degx+1):
	for j in range(degy+1):
		p = q_coeffs[i*(degy+1)+j]
		for ex in range(1,1<<f.degree()): 
			if fast_pow(beta,ex,f) == p: 
				q_coeffs_b[(i,j)] = ex
				break 
		else: 
			q_coeffs_b[(i,j)] = 0


codeword_b = []
for c in codeword: 
	for i in range(1,1<<f.degree()): 
		if c == fast_pow(beta,i,f):
			codeword_b.append(i) 
			break 
	else: 
		codeword_b.append(0) 


import sys
import subprocess

def factorize_with_sage(args):
	"""
	Calls a sage script to factorize the multivariable polynomial
	in a non-prime finite field
	"""
	a = subprocess.check_output(args)
	a = a.decode()
	return list(eval(a))

args = []
args.append("{}".format("sage"))
args.append("{}".format("factorize_with_sage.sage"))
args.append("{}".format(q_coeffs_b))
args.append("{}".format(codeword_b))
args.append("{}".format(beta.all_coeffs()[::-1]))
args.append("{}".format(f.all_coeffs()[::-1]))

msg_list_b = factorize_with_sage(args)

msg_list = []

for msg in msg_list_b:
	m = []
	for i in range(len(msg)):
		if msg[i] == 0:
			m.append(Poly(0,x,domain=FF(2)))
		else:
			m.append(fast_pow(beta,msg[i],f))
	msg_list.append(m)


print(
"""
####################
# Original Message #
####################
"""
)
print(np.array(message)[:,None])


print(
"""
#########################
# Decoded Message  List #
#########################
"""
)
for mx in msg_list:
	print(np.array(mx)[:,None])
	print()




print(
"""
#################################################################
# Q5 List Decoding Using Algorithm 2                            #
#################################################################
"""
)

n = int(input("Enter n "))
k = int(input("Enter k "))
siz = int(log(n,2))+1
f = get_irreducible(siz)
beta = find_primitive_element(f)

message = [	Poly(0,x,domain=FF(2)) for i in range(k)]

non_zero = int(input("Enter the number of non-zero coefficients "))
for i in range(non_zero):
	j = int(input("Enter power of X "))
	assert(j<len(message))
	v = int(input("Enter p, where coefficient in the form beta^p (beta is the primitive element) "))
	message[j] = fast_pow(beta,v,f)

print(
"""
####################
# Message Vector   #
####################
"""
)
print(np.array(message)[:,None])

G = gen_vandermonde(n,k,beta,f)
codeword = mat_vec_prod(G,message,f)
print(
"""
####################
# Codeword         #
####################
"""
)
print(np.array(codeword)[:,None])

dmin = n-k+1
t = (dmin-1)//2
e = int(input("Enter the number of errors (upto {}+1): ".format(t)))
assert(e <= t+1)
for i in range(e):
	j = int(input("Enter position "))
	assert(j<len(codeword))
	v = int(input("Enter p, used to derive a polynomial from F_2[x] of the form beta^p (beta is the primitive element) "))
	codeword[j] = fast_pow(beta,v,f)

print(
"""
####################
# Received Vector  #
####################
"""
)

print(np.array(codeword)[:,None])



D = ceil(sqrt(2*n*(k-1)))

constraints = []
pows = {}

for row in range(n):
	alpha = fast_pow(beta,row,f)
	y = codeword[row]
	c_row = []
	for i in range(D+1):
		for j in range(D+1):
			if i + (k-1)*j > D:
				break
			pows[(i,j)] = len(c_row)
			tmp = fast_pow(alpha,i,f) * fast_pow(y,j,f)
			_,tmp = div(tmp,f,domain=FF(2))
			c_row.append(tmp)
	c_row.append(Poly(0,x,domain=FF(2)))
	constraints.append(c_row)


q_coeffs,free = GaussElim(constraints, f, silent=True)

freei = -1
freej = -1
for i in range(D,-1,-1):
	fg = 0
	for j in range(D,-1,-1):
		if i + (k-1)*j > D:
			continue
	
		if free[pows[(i,j)]] == -1:
			freei = i
			freej = j
			fg = 1
			break
	if fg:
		break

constraints = []
for row in range(n):
	alpha = fast_pow(beta,row,f)
	y = codeword[row]
	c_row = []
	for i in range(D+1):
		for j in range(D+1):
			if i + (k-1)*j > D:
				break
			if i == freei and j == freej: 
				continue
			tmp = fast_pow(alpha,i,f) * fast_pow(y,j,f)
			_,tmp = div(tmp,f,domain=FF(2))
			c_row.append(tmp)

	tmp = fast_pow(alpha,freei,f) * fast_pow(y,freej,f)
	_,tmp = div(tmp,f,domain=FF(2))
	tmp = f-tmp
	_,tmp = div(tmp,f,domain=FF(2))
	c_row.append(tmp)
	constraints.append(c_row)

q_coeffs,_ = GaussElim(constraints, f)
q_coeffs.insert(pows[(freei,freej)], Poly(1,x,domain=FF(2)))


q_coeffs_b = {}

for i in range(D+1):
	for j in range(D+1):
		if i + (k-1)*j > D:
			break
		p = q_coeffs[pows[(i,j)]]
		for ex in range(1,1<<f.degree()): 
			if fast_pow(beta,ex,f) == p: 
				q_coeffs_b[(i,j)] = ex
				break 
		else: 
			q_coeffs_b[(i,j)] = 0


codeword_b = []
for c in codeword: 
	for i in range(1,1<<f.degree()): 
		if c == fast_pow(beta,i,f):
			codeword_b.append(i) 
			break 
	else: 
		codeword_b.append(0) 


args = []
args.append("{}".format("sage"))
args.append("{}".format("factorize_with_sage.sage"))
args.append("{}".format(q_coeffs_b))
args.append("{}".format(codeword_b))
args.append("{}".format(beta.all_coeffs()[::-1]))
args.append("{}".format(f.all_coeffs()[::-1]))

msg_list_b = factorize_with_sage(args)

msg_list = []

for msg in msg_list_b:
	m = []
	for i in range(len(msg)):
		if msg[i] == 0:
			m.append(Poly(0,x,domain=FF(2)))
		else:
			m.append(fast_pow(beta,msg[i],f))
	msg_list.append(m)


print(
"""
####################
# Original Message #
####################
"""
)
print(np.array(message)[:,None])


print(
"""
#########################
# Decoded Message  List #
#########################
"""
)
for mx in msg_list:
	print(np.array(mx)[:,None])
	print()

