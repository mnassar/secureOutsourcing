'''
Created on Apr 16, 2014

@author: NASSAR
'''
import cloud 
from paillier_gmpy2 import * 
import pickle 
import numpy as np 

path_shared="sharedWithCloud/" 
path_local="local/" 
volume_name="linearEquationsIterative"
volume_=volume_name+":"
path_cloud="/home/picloud/"+volume_name+"/"
matrixA="matrixA.npy"
vectorb="vectorb.npy"

# setup : generate the keys and send a copy of the public key to the cloud 
def setup():
    print "Generating keypair... %d bits" % 512
    priv, pub = generate_keypair(512)
    f_priv=open(path_local+"priv",'wb')
    f_pub=open(path_shared+"pub",'wb')
    f_priv.write(pickle.dumps(priv))
    f_pub.write(pickle.dumps(pub))
    f_priv.close()
    f_pub.close()  
    try: 
        cloud.volume.sync(path_shared, volume_)
        print cloud.volume.ls(volume_) 
    except cloud.cloud.CloudException:
        # executed only the first time if the volume is not already created 
        cloud.volume.create(volume_name, "/home/picloud/linearEquationsIterative", _env="with-gmpy")
        cloud.volume.sync(path_shared, volume_)
        print cloud.volume.ls(volume_)

#setup()
# load the input problem 
A= np.load(path_local+matrixA)
b= np.load(path_local+vectorb)
print "A=",A 
print "b=",b
print "find x, Ax=b"
#
#print "solve locally:"
#print np.linalg.solve(A,b)
#print "non-secure outsourcing:"
#jid=cloud.call(np.linalg.solve, A, b)
#print cloud.result(jid)
#print "secure outsourcing"

print "problem transformation to iterative form"
maxint = max(np.max(A), np.max(b))*100
size=len(b) 
D=np.diag(np.diag(A))
print "D =", D
R=A-D 
print "R=", R 
iD=np.linalg.inv(D)
c=iD.dot(b)
print "c=", c
T=-iD.dot(R)
print "T=", T
y0=np.random.randint(0, maxint, size)
y1=np.zeros(size)
eps=np.array([0.01,0.01])
print "solving the iterative form: y=Ty+c"
ite=1
while (np.abs(y1-y0)>eps).any():
    print "iter %d:" % ite,  
    y1=y0
    y0=T.dot(y0)+c
    print "delta:", y1-y0
    ite+=1
print "sol= " , y0

print "transforming and encrypting iterative form for outsourcing"
transf_factor= np.random.randint(0, maxint, size)
print "transf_factor=", transf_factor 
b1= (b+ A.dot(transf_factor))
print "b1 =",b1
c1=iD.dot(b1)
print "c1 =",c1

f_pub=open(path_shared+"pub",'rb')
pub=pickle.load(f_pub)
f_pub.close()



print "load private key for local decryption .."
f_priv=open(path_local+"priv",'rb')
priv=pickle.load(f_priv)
f_priv.close()

def fix_negative(v):
    if v-middle>0 : 
#            print "negative"
            return-(largest-v)
    else: 
        return v

print "define some vectorized functions .."
v_mpz= np.vectorize(lambda x: mpz(round(x)))
v_encrypt = np.vectorize(lambda x: encrypt(pub,mpz(round(x))))
v_decrypt = np.vectorize(lambda x: decrypt(priv, pub, x))
v_fix_negative=np.vectorize(fix_negative)



def edot(X,y):
            
    res=v_mpz(np.ones(size))
    for i in range(size):
        x=X[i]
        for j in range(size): 
            s = e_mul_const(pub, x[j], mpz(int(y[j])))
            res[i] = e_add(pub, res[i], s)
    return res

middle=pub.n/2
largest=mpz(pub.n)

scale=10**6
print "encrypt T (locally)"
T=T 
print T
eT=v_encrypt(T*scale);
print "eT=", eT
y0=np.random.randint(0, maxint, size)*scale
y1=np.zeros(size)
eps=np.array([100,100])
v_c1=v_mpz(c1*scale)
ite=1
while (np.abs(y1-y0)>eps).any():
    print "iter, %d:"% ite
    y1=y0
#    print "y1=", y1
    # equivalent of ey0=eT.dot(y0)   
    ey0=edot(eT, y0) 
    
#    print "decrytion of ey0"
    y0=v_decrypt(ey0)
    y0=v_fix_negative(y0)    
    y0/= scale
    y0+=v_c1 
#    print "y0=", y0
    
#    print "y1=", y1
#    print "y1-y0=", y1-y0
    ite+=1
print "sol= " , (y0-transf_factor*scale)/float(scale)

print "outsourcing to the cloud"
def edotCloud(y):
    f_pub=open(path_cloud+"pub",'rb')
    pub=pickle.load(f_pub)
    f_pub.close()
    X=np.load(path_cloud+"et.npy")
    res=[mpz(1)]*size
        
    for i in range(size):
        x=X[i]
        for j in range(size): 
            s = e_mul_const(pub, x[j], mpz(int(y[j])))
            res[i] = e_add(pub, res[i], s)
    return res

#np.save(path_shared+"et", eT)
#cloud.volume.sync(path_shared, volume_)
print cloud.volume.ls(volume_)

y0=np.random.randint(0, maxint, size)*scale
y1=np.zeros(size)
eps=np.array([100,100])
v_c1=v_mpz(c1*scale)
ite=1
while (np.abs(y1-y0)>eps).any():
    print "iter, %d:"% ite
    y1=y0
#    print "y1=", y1
    # equivalent of ey0=eT.dot(y0)   
    jid=cloud.call(edotCloud, y0, _vol=volume_name, _env="with-gmpy")
    ey0=cloud.result(jid)
    print ey0
#    print "decrytion of ey0"
    y0=v_decrypt(ey0)
    y0=v_fix_negative(y0)    
    y0/= scale
    y0+=v_c1 
#    print "y0=", y0
    
#    print "y1=", y1
#    print "y1-y0=", y1-y0
    ite+=1
print "sol= " , (y0-transf_factor*scale)/float(scale)




#I=np.eye(len(A))
#print "I =", I


#9.              r=randVector(b.length());    
#10.              b1=b+A*r;
#11.           D=diag(A);
#12.             R=A-D;
#13.            c=Inv(D)*b1;
#14.            T=-Inv(D)*R;
#15.            y0=randVector(b.length());   
#16.            eT=H.encrypt(T);
#17.            S.put(eT, y0); 
#18.            y1=zeros(b.length()) ;  
#19.        Run:
#20.            While(abs(y1-y0)<epsilon):
#21.                  J=S.call(H.eval, expr = eT*y0 );
#22.                  join(J) ;
#23.                y1=H.decrypt(J.result()) +c ;
#24.                S.y0=S.put(y1);
#25.            Sol=y1-r;



