import os.path #  function os.path.isfile(fname)
import random
import cPickle

class SymQuotSing(object):
	def __init__(self, r, dim=4):
		self.r = r
 		self.dim = dim
		self.ord = r**(dim-1)		#only difference here
		self.__L = ZZ**self.dim
	
		if self.dim ==2:
			self.__Q = PolynomialRing(QQ, 2, 'xy')	
		elif(self.dim == 3):
			self.__Q = PolynomialRing(QQ, 3, 'xyz')
		elif(self.dim == 4):
			self.__Q = PolynomialRing(QQ, 4, 'xyzt')
		elif(self.dim ==5):
			self.__Q = PolynomialRing(QQ, 5, 'xyztw')	
		else:
			self.__Q = PolynomialRing(QQ, self.dim, 'x')
		self.__Q.inject_variables()

	def __str__(self):
		return "Quotient of CC^"+str(self.dim)+" by the group (ZZ/"+str(self.r)+")^"+str(self.dim-1)
#different

	createDir('__EigSps')
	createDir('__ASets')
	createDir('__Relations')
	createDir('__AHilb')

	def filename_str(self):
		return "sym-"+str(self.dim)+'-'+str(self.r)

	def __ZBasis(self):
		a = []
		c = self.__L.basis()
		for i in range(self.dim):
			a.append(self.r*c[i])
		return a


#recursion to compute all the lattice points in L (I think)
        def __latptRec(self,current_vect, remaining_n, remaining_r, basket):
		if remaining_n == 0 and remaining_r == 0:
			basket.append(current_vect)
			return
		if remaining_n < 0 or remaining_r <0:
			return

	#	print current_vect, remaining_n, remaining_r, "\n"
		for i in range(remaining_r+1):
			vs = current_vect + [i]
			self.__latptRec(vs, remaining_n - 1, remaining_r -i, basket)

		return

#different
	def LatticePts(self): 	#, below = false):
		S = []
		self.__latptRec([], self.dim, self.r, S)
		
		S.sort(); 
		return S


	def __weight(self, vect):
		wt = 0
		if hasattr(vect, "exponents"):
			vect = vect.exponents()[0]
		mult = 1
		for i in range(self.dim-1,0,-1):
			wt += ((vect[0] - vect[i])%self.r)*mult
			mult *= self.r
		return wt

        def weight(self, vect):
		wt = 0
		if hasattr(vect, "exponents"):
			vect = vect.exponents()[0]
		mult = 1
		for i in range(self.dim-1,0,-1):
			wt += ((vect[0] - vect[i])%self.r)*mult
			mult *= self.r
		return wt

        

	def __eigspRecursion(self, ind, currentExponents, EigSp):
		if ind == self.dim -1:
			monomial = 1
			for k in range(self.dim):
				monomial *= self.__Q.gen(k)**currentExponents[k]

			eig = self.__weight(monomial) 

			survived = true
			for i in range(len(EigSp[eig])):
				if EigSp[eig][i] <> [0]*self.dim and greaterThan(currentExponents, EigSp[eig][i]):
					survived = false
					break
			if survived:
				EigSp[eig].append(currentExponents)
	

			return

		for j in range(self.r):
			self.__eigspRecursion(ind+1, currentExponents + [j], EigSp)	



#slightly different only in the commented row
	def EigSp(self, eigsp = [], exponents_only=false):
		
		if eigsp != []:
			return eigsp

		if os.path.exists('__EigSps/'+self.filename_str()+'__eigsps.p'):
			f = open('__EigSps/'+self.filename_str()+'__eigsps.p', 'rb')
			EigSp = cPickle.load(f)
			f.close()
			
		else:
			EigSp = []
			for i in range(self.ord):
				EigSp.append([])

			for j in range(self.dim):
				#if gcd(self.r, self.P[j])==1:
				t = [0]*self.dim
				t[j] = self.r
				EigSp[0].append(t)

			self.__eigspRecursion(-1, [], EigSp) 

			f = open('__EigSps/'+self.filename_str()+'__eigsps.p', 'wb')
			cPickle.dump(EigSp, f)
			f.close()		


		if exponents_only:
			return EigSp

		for i in range(self.ord):
			for j in range(len(EigSp[i])):
				monom = 1
				for k in range(self.dim):
					monom *= self.__Q.gen(k)**EigSp[i][j][k]
				EigSp[i][j] = monom
		return EigSp



	def ASets(self, eigsp = [], exponents_only = false):

		if os.path.exists('__ASets/'+self.filename_str()+'__a-sets.p'):
			f = open('__ASets/'+self.filename_str()+'__a-sets.p', 'rb')
			asets = cPickle.load(f)
			f.close()
			
			for i in range(len(asets)):
				for j in range(len(asets[i])):
					monom = 1
					for k in range(self.dim):
						monom *= self.__Q.gen(k)**asets[i][j][k]
					asets[i][j] = monom
			return asets


		EigSp = self.EigSp(eigsp)
		I = []
		C = []
		M = []

		asets = []
		finished = false
		while not finished:
			S = exclude(EigSp[0],1)

			over = true
			max_i = -1
			for i in range(len(I)):
				S += exclude(EigSp[I[i]], C[i][M[i]])
				if len(C[i]) != M[i]+1:
					over = false
					max_i = i

			Id = self.__Q.ideal(S)
			Qbar = self.__Q.quotient_ring(Id)

			remaining = []
			exists_empty = false
			for i in range(1, self.ord):
				surviving_i = []
				for j in range(len(EigSp[i])):
					el = EigSp[i][j]
					if el not in Id:
						surviving_i.append(el_quot(el, Qbar))
				if len(surviving_i) != 1:
					remaining.append([i,surviving_i])
				if len(surviving_i) == 0:
					exists_empty = true


			if len(remaining) == 0:
				asets.append(Qbar.defining_ideal().interreduced_basis())

			if len(remaining) != 0 and exists_empty == false:
				I.append(remaining[0][0])
				C.append(remaining[0][1])
				M.append(0)

			if len(remaining) == 0 or (len(remaining)!=0 and exists_empty ==true):
				if over:
					finished = true

				else:
					broj = len(I) - max_i - 1
					remove_last(I, broj)
					remove_last(M, broj)
					remove_last(C, broj)
					M[max_i] += 1 

		asets2 = []
		for i in range(len(asets)):
			asets2.append([])
			for j in range(len(asets[i])):
				asets2[i].append( asets[i][j].exponents()[0])

		f = open('__ASets/'+self.filename_str()+'__a-sets.p', 'wb')
		cPickle.dump(asets2, f)
		f.close()		
		
		if exponents_only:
			return asets2

		return asets

	def __eigenbasis(self, S):
		m = 1
		for i in range(self.dim):
			m *= self.__Q.gen(i)**self.r
		basis = Set( self.__Q.monomial_all_divisors(m) )

		for I in range(len(S)):
			opp = self.__Q.monomial_quotient(m, S[I])
			L = self.__Q.monomial_all_divisors(opp)
			basis -= Set(L)

	#	print "basis is  ", basis

		if self.ord != len(basis):
			print "Not a cluster!"
			return []
		
		clus = [0]*self.ord

		for i in range(self.ord):
			mono = self.__Q.monomial_quotient(m, basis[i])
			ind = self.__weight(mono)
			clus[ind] = mono
		return clus

	def EigenBases(self):

#		if os.path.exists('__EigenBases/'+self.filename_str()+'__eigbas.p'):
#			f = open('__EigenBases/'+self.filename_str()+'__eigbas.p', 'rb')
#			basket  = cPickle.load(f)
#			f.close()
#		else:

		basket = []
		A = self.ASets()
		for i in range(len(A)):
			basket.append( self.__eigenbasis(A[i]) )
#		f = open('__EigenBases/'+self.filename_str()+'__eigbas.p', 'wb')
#		cPickle.dump(basket, f)
#		f.close()	
	
		return basket



	def __reduce_rels(self, mat):
	        M = matrix(mat)
		K = M.kernel().matrix().rows()
		indices = []

		for i in range(len(K)):
			npos = 0
			nneg = 0
			indp = indn = -1
			for j in range(len(mat)):
				if K[i][j] > 0:
					npos += 1
					indp = j
				else:
					if K[i][j] < 0:
						nneg += 1
						indn = j
	#                       print i, j, npos, indp, nneg, indn
			if npos == 1  and K[i][indp] ==1 and nneg > 0:
				indices.append(indp)
			else:
				if npos > 0 and nneg == 1 and K[i][indn]==-1:
					indices.append(indn)
			#print indices
	#	print "indexes to be removed: ", indices
		M = M.delete_rows(indices)
		return M.rows()


	def __relations(self, aset, basis): #list of relations
		mat = []
		for i in range(len(aset)):
			bb = aset[i].exponents()[0]
			ind = self.__weight(bb)
			cc = basis[ind].exponents()[0]
			#print bb, cc
			row = []
			for i in range(len(bb)):
				row.append(bb[i] - cc[i])
			mat.append( row )
#		print mat
		redmat = self.__reduce_rels(mat)
		if len(redmat) != self.dim:
			random.shuffle(redmat)
		while (redmat != mat):
			mat = redmat
			redmat = self.__reduce_rels(mat)
		for i in range(self.r):
			random.shuffle(redmat)
			mat = redmat
			redmat = self.__reduce_rels(mat)
		return mat
	
	def Relations(self):

		if os.path.exists('__Relations/'+self.filename_str()+'__rels.p'):
			f = open('__Relations/'+self.filename_str()+'__rels.p', 'rb')
			basket  = cPickle.load(f)
			f.close()

		else:
			basket = []
			A = self.ASets()
			B = self.EigenBases()
			for i in range(len(A)):
				basket.append( self.__relations(A[i], B[i]) )
			f = open('__Relations/'+self.filename_str()+'__rels.p', 'wb')
			cPickle.dump(basket, f)
			f.close()	

		return basket


	def __affinepiece(self, mat):
		pts = []
		if len(mat) == self.dim:
			A = (1/self.r**(self.dim-2))*matrix(mat).adjoint()
			for i in range(self.dim):
				if A[0][i] < 0:
					A *= -1
					break
				if A[0][i] > 0:
					break
			pts = A.columns()
		return pts


	def AHilb(self):

		if os.path.exists('__AHilb/'+self.filename_str()+'__AHilb.p'):
			f = open('__AHilb/'+self.filename_str()+'__AHilb.p', 'rb')
			pieces  = cPickle.load(f)
			f.close()

		else:
			matrices = self.Relations()
			pieces = []
			for i in range(len(matrices)):
				pieces.append( self.__affinepiece(matrices[i]) )
			f = open('__AHilb/'+self.filename_str()+'__AHilb.p', 'wb')
			cPickle.dump(pieces, f)
			f.close()

		return pieces


	def __project_pt(self, P):
		G = [0]*self.dim
		for i in range(self.dim):
			G[i] = P[i]/ptAge(P)
		return G

