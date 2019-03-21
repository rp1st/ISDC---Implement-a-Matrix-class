import math
from math import sqrt
import numbers

def zeroes(height, width):         # Global function #
        """
        Creates a matrix of zeroes.
        """
        if height == 1:
            g = [[0.0]]
            return Matrix(g)
        else:
            g = [[0.0 for _ in range(width)] for __ in range(height)]
            return Matrix(g)

def identity(n):                 # Global function #
        """
        Creates a n x n identity matrix.
        """
        I = zeroes(n, n)
        for i in range(n):
            I.g[i][i] = 1.0
        return I
    
def multiply(A, B):              # Global function #
        """ 
        Multiply the matrix by a matrix or a matrix by a scalar value.
        This is a private function called by __mul__ and __rmul__.
        """
        assert isinstance(A, (Matrix, int, float)) and isinstance(B, (Matrix, int, float))
        if isinstance(A, Matrix) and isinstance(B, Matrix):
            if len(A.g[0]) != len(B.g):
                raise TypeError("Both matrices should be either a compatible Matrix object, or a constant")
            res = zeroes(len(A.g), len(B.g[0]))
            for row in range(len(A.g)):
                for col in range(len(B.g[row])):
                    for k in range(len(B.g)):
                        res.g[row][col] += A.g[row][k] * B.g[k][col]        # Dot product between row and column elements #
            return res
        left, right = (A, B) if isinstance(A, Matrix) else (B, A)
        res = [[left[col][row] * right for row in range(len(B.g[0]))]for col in range(len(B.g))]
        return Matrix(res)

    
class Matrix(object):         # Class begins from here #

    # Constructor
    def __init__(self, grid):
        self.g = grid
        self.h = len(grid)
        self.w = len(grid[0])

    #
    # Primary matrix math methods
    #############################
## 3x3 and 4x4 Determinant and Inverse formulae from http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html ##

    def determinant(self):
        """
        Calculates the determinant of a 1x1 or 2x2 matrix.
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate determinant of non-square matrix.")
        if self.h > 4:
            raise(NotImplementedError, "Calculating determinant not implemented for matrices larger than 4x4")
        
        # TODO - your code here
        # DONE - see my code below
        if self.h == 1:
            return self.g[0][0] # Determinant is the element itself #

        if self.h == 2: #Formula for 2x2 matrix determinant is a*d - b*c #
            return (self.g[0][0]*self.g[1][1])-(self.g[0][1]*self.g[1][0])
        
        if self.h == 3: #Formula for 3x3 matrix determinant is detA=a11a22a33 + a21a32a13 + a31a12a23 - a11a32a23 - a31a22a13 - a21a12a33 #
            return (self.g[0][0]*self.g[1][1]*self.g[2][2])+(self.g[1][0]*self.g[2][1]*self.g[0][2])+(self.g[2][0]*self.g[0][1]*self.g[1][2])-(self.g[0][0]*self.g[2][1]*self.g[1][2])-(self.g[2][0]*self.g[1][1]*self.g[0][2])-(self.g[1][0]*self.g[0][1]*self.g[2][2])
            
        if self.h == 4: # Formula for 4x4 is big and will not fit in one line #
            return (self.g[0][0]*self.g[1][1]*self.g[2][2]*self.g[3][3])+(self.g[0][0]*self.g[1][2]*self.g[2][3]*self.g[3][1])+(self.g[0][0]*self.g[1][3]*self.g[2][1]*self.g[3][2])+(self.g[0][1]*self.g[1][0]*self.g[2][3]*self.g[3][2])+(self.g[0][1]*self.g[1][2]*self.g[2][0]*self.g[3][3])+(self.g[0][1]*self.g[1][3]*self.g[2][2]*self.g[3][0])+(self.g[0][2]*self.g[1][0]*self.g[2][1]*self.g[3][3])+(self.g[0][2]*self.g[1][1]*self.g[2][3]*self.g[3][0])+(self.g[0][2]*self.g[1][3]*self.g[2][0]*self.g[3][1])+(self.g[0][3]*self.g[1][0]*self.g[2][2]*self.g[3][1])+(self.g[0][3]*self.g[1][1]*self.g[2][0]*self.g[3][2])+(self.g[0][3]*self.g[1][2]*self.g[2][1]*self.g[3][0])-(self.g[0][0]*self.g[1][1]*self.g[2][3]*self.g[3][2])-(self.g[0][0]*self.g[1][2]*self.g[2][1]*self.g[3][3])-(self.g[0][0]*self.g[1][3]*self.g[2][2]*self.g[3][1])-(self.g[0][1]*self.g[1][0]*self.g[2][2]*self.g[3][3])-(self.g[0][1]*self.g[1][2]*self.g[2][3]*self.g[3][0])-(self.g[0][1]*self.g[1][3]*self.g[2][0]*self.g[3][2])-(self.g[0][2]*self.g[1][0]*self.g[2][3]*self.g[3][1])-(self.g[0][2]*self.g[1][1]*self.g[2][0]*self.g[3][3])-(self.g[0][2]*self.g[1][3]*self.g[2][1]*self.g[3][0])-(self.g[0][3]*self.g[1][0]*self.g[2][1]*self.g[3][2])-(self.g[0][3]*self.g[1][1]*self.g[2][2]*self.g[3][0])-(self.g[0][3]*self.g[1][2]*self.g[2][0]*self.g[3][1])
            
            
    def trace(self):
        """
        Calculates the trace of a matrix (sum of diagonal entries).
        """
        if not self.is_square():
            raise(ValueError, "Cannot calculate the trace of a non-square matrix.")

        # TODO - your code here
        # DONE - see my code below
        tr = 0
        for i in range(self.h):
            tr += self.g[i][i]        # Formula for computing trace #
        return tr
    
    
    def inverse(self):
        """
        Calculates the inverse of a 1x1, 2x2, 3x3 and 4x4 Matrix.
        """
        if not self.is_square():
            raise(ValueError, "Non-square Matrix does not have an inverse.")
        if self.h > 4:
            raise(NotImplementedError, "inversion not implemented for matrices larger than 4x4")

        # TODO - your code here
        # DONE - see my code below
        if self.h == 1:        # for 1x1 matrix #
            inverse = zeroes(self.h, self.g)
            inverse.g[0][0] = 1/self.g[0][0]       # Formula for 1x1 matrix #
            return inverse

        elif self.h == 2:      # for 2x2 matrix #
            inverse = zeroes(self.h, self.w)
            det = self.determinant()
            if det != 0:
                a = 1/det
                inverse[0][0] = a*self.g[1][1]          # Formulae for 2x2 matrix #
                inverse[0][1] = a*(-1)*self.g[0][1]     #
                inverse[1][0] = a*(-1)*self.g[1][0]     #
                inverse[1][1] = a*self.g[0][0]          #
                return inverse
            else:
                raise ValueError("No inverse, since the determinant of matrix is zero")
                return 
            
        elif self.h == 3:       # for 3x3 matrix #
            inverse = zeroes(self.h, self.w)
            det = self.determinant()
            if det != 0:
                a = 1/det
                inverse[0][0] = a*(self.g[1][1]*self.g[2][2] - self.g[1][2]*self.g[2][1])     # Formulae for 3x3 matrix #
                inverse[0][1] = a*(self.g[0][2]*self.g[2][1] - self.g[0][1]*self.g[2][2])     #
                inverse[0][2] = a*(self.g[0][1]*self.g[1][2] - self.g[0][2]*self.g[1][1])     #
                inverse[1][0] = a*(self.g[1][2]*self.g[2][0] - self.g[1][0]*self.g[2][2])     #
                inverse[1][1] = a*(self.g[0][0]*self.g[2][2] - self.g[0][2]*self.g[2][0])     #
                inverse[1][2] = a*(self.g[0][2]*self.g[1][0] - self.g[0][0]*self.g[1][2])     #
                inverse[2][0] = a*(self.g[1][0]*self.g[2][1] - self.g[1][1]*self.g[2][0])     #
                inverse[2][1] = a*(self.g[0][1]*self.g[2][0] - self.g[0][0]*self.g[2][1])     #
                inverse[2][2] = a*(self.g[0][0]*self.g[1][1] - self.g[0][1]*self.g[1][0])     #
                return inverse
            else:
                raise ValueError("No inverse, since the determinant of matrix is zero")
                return
            
        elif self.h == 4:       # for 4x4 matrix #
            inverse = zeroes(self.h, self.w)
            det = self.determinant()
            if det != 0:
                a = 1/det                   # Formulae for 4x4 matrix #
                inverse[0][0] = a*(self.g[1][1]*self.g[2][2]*self.g[3][3]+self.g[1][2]*self.g[2][3]*self.g[3][1]+self.g[1][3]*self.g[2][1]*self.g[3][2]-self.g[1][1]*self.g[2][3]*self.g[3][2]-self.g[1][2]*self.g[2][1]*self.g[3][3]-self.g[1][3]*self.g[2][2]*self.g[3][1])
                inverse[0][1] = a*(self.g[0][1]*self.g[2][3]*self.g[3][2]+self.g[0][2]*self.g[2][1]*self.g[3][3]+self.g[0][3]*self.g[2][2]*self.g[3][1]-self.g[0][1]*self.g[2][2]*self.g[3][3]-self.g[0][2]*self.g[2][3]*self.g[3][1]-self.g[0][3]*self.g[2][1]*self.g[3][2])
                inverse[0][2] = a*(self.g[0][1]*self.g[1][2]*self.g[3][3]+self.g[0][2]*self.g[1][3]*self.g[3][1]+self.g[0][3]*self.g[1][1]*self.g[3][2]-self.g[0][1]*self.g[1][3]*self.g[3][2]-self.g[0][2]*self.g[1][1]*self.g[3][3]-self.g[0][3]*self.g[1][2]*self.g[3][1])
                inverse[0][3] = a*(self.g[0][1]*self.g[1][3]*self.g[2][2]+self.g[0][2]*self.g[1][1]*self.g[2][3]+self.g[0][3]*self.g[1][2]*self.g[2][1]-self.g[0][1]*self.g[1][2]*self.g[2][3]-self.g[0][2]*self.g[1][3]*self.g[2][1]-self.g[0][3]*self.g[1][1]*self.g[2][2])                     
                inverse[1][0] = a*(self.g[1][0]*self.g[2][3]*self.g[3][2]+self.g[1][2]*self.g[2][0]*self.g[3][3]+self.g[1][3]*self.g[2][2]*self.g[3][0]-self.g[1][0]*self.g[2][2]*self.g[3][3]-self.g[1][2]*self.g[2][3]*self.g[3][0]-self.g[1][3]*self.g[2][0]*self.g[3][2])
                inverse[1][1] = a*(self.g[0][0]*self.g[2][2]*self.g[3][3]+self.g[0][2]*self.g[2][3]*self.g[3][0]+self.g[0][3]*self.g[2][0]*self.g[3][2]-self.g[0][0]*self.g[2][3]*self.g[3][2]-self.g[0][2]*self.g[2][0]*self.g[3][3]-self.g[0][3]*self.g[2][2]*self.g[3][0])
                inverse[1][2] = a*(self.g[0][0]*self.g[1][3]*self.g[3][2]+self.g[0][2]*self.g[1][0]*self.g[3][3]+self.g[0][3]*self.g[1][2]*self.g[3][0]-self.g[0][0]*self.g[1][2]*self.g[3][3]-self.g[0][2]*self.g[1][3]*self.g[3][0]-self.g[0][3]*self.g[1][0]*self.g[3][2])
                inverse[1][3] = a*(self.g[0][0]*self.g[1][2]*self.g[2][3]+self.g[0][2]*self.g[1][3]*self.g[2][0]+self.g[0][3]*self.g[1][0]*self.g[2][2]-self.g[0][0]*self.g[1][3]*self.g[2][2]-self.g[0][2]*self.g[1][0]*self.g[2][3]-self.g[0][3]*self.g[1][2]*self.g[2][0])                     
                inverse[2][0] = a*(self.g[1][0]*self.g[2][1]*self.g[3][3]+self.g[1][1]*self.g[2][3]*self.g[3][0]+self.g[1][3]*self.g[2][0]*self.g[3][1]-self.g[1][0]*self.g[2][3]*self.g[3][1]-self.g[1][1]*self.g[2][0]*self.g[3][3]-self.g[1][3]*self.g[2][1]*self.g[3][0])
                inverse[2][1] = a*(self.g[0][0]*self.g[2][3]*self.g[3][1]+self.g[0][1]*self.g[2][0]*self.g[3][3]+self.g[0][3]*self.g[2][1]*self.g[3][0]-self.g[0][0]*self.g[2][1]*self.g[3][3]-self.g[0][1]*self.g[2][3]*self.g[3][0]-self.g[0][3]*self.g[2][0]*self.g[3][1])
                inverse[2][2] = a*(self.g[0][0]*self.g[1][1]*self.g[3][3]+self.g[0][1]*self.g[1][3]*self.g[3][0]+self.g[0][3]*self.g[1][0]*self.g[3][1]-self.g[0][0]*self.g[1][3]*self.g[3][1]-self.g[0][1]*self.g[1][0]*self.g[3][3]-self.g[0][3]*self.g[1][1]*self.g[3][0])
                inverse[2][3] = a*(self.g[0][0]*self.g[1][3]*self.g[2][1]+self.g[0][1]*self.g[1][0]*self.g[2][3]+self.g[0][3]*self.g[1][1]*self.g[2][0]-self.g[0][0]*self.g[1][1]*self.g[2][3]-self.g[0][1]*self.g[1][3]*self.g[2][0]-self.g[0][3]*self.g[1][0]*self.g[2][1])                     
                inverse[3][0] = a*(self.g[1][0]*self.g[2][2]*self.g[3][1]+self.g[1][1]*self.g[2][0]*self.g[3][2]+self.g[1][2]*self.g[2][1]*self.g[3][0]-self.g[1][0]*self.g[2][1]*self.g[3][2]-self.g[1][1]*self.g[2][2]*self.g[3][0]-self.g[1][2]*self.g[2][0]*self.g[3][1])
                inverse[3][1] = a*(self.g[0][0]*self.g[2][1]*self.g[3][2]+self.g[0][1]*self.g[2][2]*self.g[3][0]+self.g[0][2]*self.g[2][0]*self.g[3][1]-self.g[0][0]*self.g[2][2]*self.g[3][1]-self.g[0][1]*self.g[2][0]*self.g[3][2]-self.g[0][2]*self.g[2][1]*self.g[3][0])
                inverse[3][2] = a*(self.g[0][0]*self.g[1][2]*self.g[3][1]+self.g[0][1]*self.g[1][0]*self.g[3][2]+self.g[0][2]*self.g[1][1]*self.g[3][0]-self.g[0][0]*self.g[1][1]*self.g[3][2]-self.g[0][1]*self.g[1][2]*self.g[3][0]-self.g[0][2]*self.g[1][0]*self.g[3][1])
                inverse[3][3] = a*(self.g[0][0]*self.g[1][1]*self.g[2][2]+self.g[0][1]*self.g[1][2]*self.g[2][0]+self.g[0][2]*self.g[1][0]*self.g[2][1]-self.g[0][0]*self.g[1][2]*self.g[2][1]-self.g[0][1]*self.g[1][0]*self.g[2][2]-self.g[0][2]*self.g[1][1]*self.g[2][0])
                return inverse
            else:
                raise ValueError("No inverse, since the determinant of matrix is zero")
                return
        else:
            return 

 
    def T(self):
        """
        Returns a transposed copy of this Matrix.
        """
        # TODO - your code here
        # DONE - see my code below
        return Matrix([[row[i] for row in self] for i in range(self.w)])   
    
    #
    # Begin Operator Overloading
    ############################
    def __getitem__(self,idx):
        """
        Defines the behavior of using square brackets [] on instances
        of this class.

        Example:

        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > my_matrix[0]
          [1, 2]

        > my_matrix[0][0]
          1
        """
        return self.g[idx]

    
    def __repr__(self):
        """
        Defines the behavior of calling print on an instance of this class.
        """
        s = ""
        for row in self.g:
            s += " ".join(["{} ".format(x) for x in row])
            s += "\n"
        return s

    
    def __add__(self,other):
        """
        Defines the behavior of the + operator
        """
        if self.h != other.h or self.w != other.w:
            raise(ValueError, "Matrices can only be added if the dimensions are the same") 
        #   
        # TODO - your code here
        # DONE - see my code below
        res = zeroes(self.h, self.w)         # Initialize zero resultant matrix #
        for row in range(self.h):
            for col in range(self.w):
                res.g[row][col] = self.g[row][col] + other.g[row][col]         # Formula to add individual elements of 2 matrices #
        return res

    def __neg__(self):
        """
        Defines the behavior of - operator (NOT subtraction)
        Example:
        > my_matrix = Matrix([ [1, 2], [3, 4] ])
        > negative  = -my_matrix
        > print(negative)
          -1.0  -2.0
          -3.0  -4.0
        """
        # TODO - your code here
        # DONE - see my code below
        res = zeroes(self.h, self.w)          # Initialize zero resultant matrix #
        for row in range(self.h):
            for col in range(self.w):
                res.g[row][col] = self.g[row][col] * (-1)         # Multiply individual elements of matrix with -1 #
        return res
    
        
    def __sub__(self, other):
        """
        Defines the behavior of - operator (as subtraction)
        """
        #   
        # TODO - your code here
        # DONE - see my code below
        res = zeroes(self.h, self.w)           # Initialize zero resultant matrix #
        for row in range(self.h):
            for col in range(self.w):
                res.g[row][col] = self.g[row][col] - other.g[row][col]    # Formula to substract individual elements of 2 matrices #
        return res


    def __mul__(self, other):
        """
        Defines the behavior of * operator (matrix multiplication)
        """
        #   
        # TODO - your code here
        # DONE - see my code below
        return multiply(self, other)         # Calls global function multiply() for multiplication of two matrices #

    
    def __rmul__(self, other):
        """
        Called when the thing on the left of the * is not a matrix.
        Example:
        > identity = Matrix([ [1,0], [0,1] ])
        > doubled  = 2 * identity
        > print(doubled)
          2.0  0.0
          0.0  2.0
        """
        if isinstance(other, numbers.Number):
            pass
            #   
            # TODO - your code here
            # DONE - see my code below
            return multiply(other, self)        # Calls global function multiply() for multiplication of scalar with matrix #

        
    def is_square(self):
        return self.h == self.w
