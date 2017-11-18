def createDir(filename):
	try:
		if not os.path.exists(filename):
			os.makedirs(filename)
	except OSError:
		print "cannot"

#note - strictly greater
def greaterThan(vector1, vector2):
        n = len(vector1)
        try:
                for i in range(n):
                        if (vector1[i] < vector2[i]):
                                return false
                return true
        except IndexError:
                print("Error: vectors of different dimensions")


def ptAge(point):
        sum = 0
        for i in range(len(point)):
                sum += point[i]
        return sum

def contFr(r,a):
        c = []
        while  a > 0 :
                c_i = ceil(r/a)
                c.append(c_i)
                r_old = r
                r = a
                a = c_i * a - r_old
        return c

def frCont(S):
        n = len(S)
        if (n==1): return S[0]

        denom = 1
        num = S[n-1]
        i = n-1
        while(i > 0):
                temp = denom
                denom = num
                num = S[i-1]*num - temp
                i-=1
        return num, denom


def ModInverse(s,r):
        if gcd(s,r) > 1:
                return 0
        return inverse_mod(s,r)


#returns the list 'lis' without the element 'element' (does not change the original list)
def exclude(lis, element):
	a = []
	for i in range(len(lis)):
		if lis[i] != element:
			a.append(lis[i])
	return a


#removes last n elements from the list
def remove_last(lis, n):
	for i in range(n):
		lis.pop()


def el_quot(m, A):
	return A.lift(A.retract(m))

def decToBase(num,base):
    a =0
    i  =1
    while(num):
	rem = num % base
	a = a+ i*rem
	num = (num - rem) /  base
	i = i*10
    return a

