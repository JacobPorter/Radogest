
class A:
    def __init__(self, i):
        self.i = i

class B(A):
    def b_add(self):
        self.i += 1

class C(A):
    def c_add(self):
        self.i += 1

class D(B, C):
    def test(self):
        self.c_add()
        print(self.i)
        self.b_add()
        print(self.i)

my_D = D(0)
print(isinstance(my_D, A))
print(my_D.i)
my_D.test()
my_D.test()
