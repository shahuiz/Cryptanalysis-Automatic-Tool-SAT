import numpy as np
import subprocess
import time

# Cipher Details
BlockSize = 64
FullRound = 31
SBoxLength = 4
SBoxSize = 16
NumOfParallelSBoxes = 16
S_Box = [12, 5, 6, 11, 9, 0, 10, 13, 3, 14, 15, 8, 4, 7, 1, 2]
P_Box = [0, 16, 32, 48, 1, 17, 33, 49, 2, 18, 34, 50, 3, 19, 35, 51, 4, 20, 36, 52, 5, 21, 37, 53, 6, 22, 38, 54, 7, 23, 39, 55, 8, 24, 40, 56, 9, 25, 41, 57, 10, 26, 42, 58, 11, 27, 43, 59, 12, 28, 44, 60, 13, 29, 45, 61, 14, 30, 46, 62, 15, 31, 47, 63]

# Generate DDT for PRESENT S-Box
ddt = np.zeros((16, 16))
truncated_ddt = np.zeros((16, 16))
for input_diff in range (0, 16):
   for x in range (0, 16):
      output_diff = S_Box[x] ^ S_Box[x ^ input_diff]
      ddt[input_diff][output_diff] += 1
      truncated_ddt[input_diff][output_diff] = 1

# !!IMPORTANT!! Espresso I/O is in SOP, while SAT solver forces input in POS
# Inverse DDT: 0-> possible 1-> impossible, for later conversion from SOP to POS
inv_trunc_ddt = np.int64(truncated_ddt == 0)

# Generate PRESENT S-Box Constraints in SOP, 1 for impossible pattern, 0 for possible pattern
pladoc = []
pladoc.append('.i ' + str(2 * SBoxLength + 1))
pladoc.append('.o 1')
for i in range(16):
   for j in range(16):            
      diff = str(bin(i)[2:].zfill(SBoxLength)) + str(bin(j)[2:].zfill(SBoxLength))
      if(diff[0:4] == '0000'):   # Inactive S-Box
         pladoc.append(diff + "0 " + str(inv_trunc_ddt[i][j]))
         pladoc.append(diff + "1 1")
      else:                      # Active S-Box
         pladoc.append(diff + "0 1")
         pladoc.append(diff + "1 " + str(inv_trunc_ddt[i][j]))
pladoc.append(".e")

# Feed the PLA file to ESPRESSO to minimize clause count
result = subprocess.run(['espresso'], input = '\n'.join(pladoc), capture_output=True, text=True)

# Crop and parse the result file to S-Box Constraints
SBoxConstraints = result.stdout.split("\n")[3:-2]
# Iterpretation of ESPRESSO result: (SOP) '1' = ON, '0' = OFF, '-' = DONT_CARE; only return result = 1, for impossible patterns

# Convert SOP form to POS form: ON(SOP)->OFF(POS), OFF(SOP)->ON(POS), DONT_CARE(SOP)->DONT_CARE(POS)
# Effect: S-Box constraints are each generated to exclude impossible patterns (result = 1 in SOP)
def GenSBoxConstraint(cnfdoc, x):
   for i in range(len(SBoxConstraints)):
      tempclause = ''
      for j in range(9):
         if(SBoxConstraints[i][j] == '-'):
            continue
         if(SBoxConstraints[i][j] == '1'):
            tempclause += '-'
         tempclause = tempclause + str(x[j]) + ' '
      cnfdoc.append(tempclause + '0')

def GenObjectiveFunction(cnfdoc, x, s, n, k):
   # Default
   if (k == 0):
      for i in range(n):
         cnfdoc.append('-' + str(x[i]) + ' 0')
      return
   # Generate Sequential Counter Circuit
   cnfdoc.append('-' + str(x[0]) + ' ' + str(s[0][0]) + ' 0')
   for j in range(1, k):
      cnfdoc.append('-' + str(s[0][j]) + ' 0')
   for i in range(1, n-1):
      cnfdoc.append('-' + str(x[i]) + ' ' + str(s[i][0]) + ' 0')
      cnfdoc.append('-' + str(s[i-1][0]) + ' ' + str(s[i][0]) + ' 0')
   for j in range(1, k):
      for i in range(1, n-1):
         cnfdoc.append('-' + str(x[i]) + ' ' + '-' + str(s[i-1][j-1]) + ' ' + str(s[i][j]) + ' 0')
         cnfdoc.append('-' + str(s[i-1][j]) + ' ' + str(s[i][j]) + ' 0')
   for i in range(1, n-1):
      cnfdoc.append('-' + str(x[i]) + ' ' + '-' + str(s[i-1][k-1]) + ' 0')
   cnfdoc.append('-' + str(x[n-1]) + ' ' + '-' + str(s[n-2][k-1]) + ' 0')
   return

def GenMatsuiSetOfBound(type, targetround, totalround):
   MatsuiBounds = []
   if (type == "S"):
      start = targetround
      for end in range(start + 1, totalround + 1):
         MatsuiBounds.append([start, end])
   if (type == "E"):
      end = targetround
      if (end > totalround):
         end = totalround
      for start in range(end-1, -1, -1):
         MatsuiBounds.append([start, end])
   return MatsuiBounds

def GenMatsuiConstraint(cnfdoc, x, s, n, k, left, right, m):
   # Default
   if (m == 0):
      for i in range(left, right + 1):
         cnfdoc.append('-' + str(x[i]) + ' 0')
      return
   # Generate Matsui's Bounding Condition 
   if (left == 0 and right < n-1):
      for i in range(1, right + 1):
         cnfdoc.append('-' + str(x[i]) + ' ' + '-' + str(s[i-1][m-1]) + ' 0')
   if (left > 0 and right == n-1):
      for i in range(0, k-m):
         cnfdoc.append(str(s[left-1][i]) + ' ' + '-' + str(s[n-2][i+m]) + ' 0')
      for i in range(0, k-m+1):
         cnfdoc.append(str(s[left-1][i]) + ' ' + '-' + str(x[n-1]) + ' ' + '-' + str(s[n-2][i+m-1]) + ' 0')
   if (left > 0 and right < n-1):
      for i in range(0, k-m):
         cnfdoc.append(str(s[left-1][i]) + ' ' + '-' + str(s[right][i+m]) + ' 0')
   return

def FormulateCNF(round, k, MatsuiBounds):
   # Allocate varaiables
   NumOfStateVar = (round + 1) * BlockSize
   differential = np.arange(1, NumOfStateVar + 1).reshape(round + 1, BlockSize)
   n = round * NumOfParallelSBoxes
   roundin = differential[0: round]
   roundout = differential[1 : round + 1]   
   # w: dummy variables record the activeness of S-Boxes at each round, objective variables
   NumOfDummyVar = round * NumOfParallelSBoxes 
   w = np.arange(NumOfStateVar+ 1, NumOfStateVar + NumOfDummyVar + 1).reshape(round, NumOfParallelSBoxes)
   # s: auxiliary variables to be used in sequential counter circuit, record partial sums
   NumOfSCCVar = (NumOfDummyVar - 1) * k
   s = np.arange(NumOfStateVar + NumOfDummyVar + 1, NumOfStateVar + NumOfDummyVar + NumOfSCCVar + 1).reshape(NumOfDummyVar -1, k)
   # Initialize DIMAC file
   cnfdoc = ['p cnf ', '']
   # Force nonzero input
   for i in range(BlockSize):
      cnfdoc[1] += str(roundin[0][i]) + ' '
   cnfdoc[1] += '0'
   # Loop all rounds
   for r in range(round):
      # Denote round function as input = a -S-Box-> b -P-Box-> c = output
      a = roundin[r]
      c = roundout[r]
      b = np.zeros(BlockSize, np.int64)
      for i in range(BlockSize):
         b[i] = c[P_Box[i]]
      # Enumerate vector x = a||b||w for each S-Box
      for i in range(NumOfParallelSBoxes):
         x = np.concatenate((a[4*i : 4*(i+1)], b[4*i : 4*(i+1)], np.asarray(w[r][i])), axis = None)
         GenSBoxConstraint(cnfdoc, x)
      # Generate objective function
      objvar = w.reshape(-1)
      GenObjectiveFunction(cnfdoc, objvar, s, n, k)
      # Generate Matsui's Bounding Conditions
      for [startround, endround] in MatsuiBounds:
         left, right = 16 * startround, 16 * endround - 1
         m = k - DiffActiveSbox[startround] - DiffActiveSbox[round - endround]
         GenMatsuiConstraint(cnfdoc, objvar, s, n, k, left, right, m)
   # Complete DIMAC file
   TotalNumOfVar = NumOfStateVar + NumOfDummyVar + NumOfSCCVar 
   TotalNumOfClause = len(cnfdoc) - 1
   cnfdoc[0] += str(TotalNumOfVar) + ' ' + str(TotalNumOfClause)
   return cnfdoc

def CallSolver(cnfdoc, solver):
   if (solver == "Cryptominisat"):
      order = ["Cryptominisat5", "--verb", "0"]
   if (solver == "Cadical"):
      order = ["cadical", "--sat", "-q"]
   
   start_time = time.time()
   solution = subprocess.run(order, input ='\n'.join(cnfdoc), capture_output=True, text=True).stdout
   if(solution.startswith("s SATISFIABLE")):
      return time.time() - start_time, solution
   else:
      return -1, ""

# Start of main function
# Initialization
DiffActiveSbox = np.zeros(FullRound + 1, np.int64)
SolutionSet = []

# Ask for input
solver_choice = str(input("Enter choice of solver, Cryptominisat/Cadical:> "))
type_matsui = str(input("Enter type of Matsui's bounding condition, S for StartAt and E for EndAt:> "))
target_round = int(input("Enter target round for Matsui's bounding condition, integer between 0 and 31:> "))

# Clear run log
open("./runlog_" + solver_choice + '_' + type_matsui + '_' + str(target_round), "w").close()

total_start = time.time()
# Search start
for round in range(1, 32):
   MatsuiBounds =  GenMatsuiSetOfBound(type_matsui, target_round, round)
   SBoxLimit = DiffActiveSbox[round-1] + 1   # Assmption: After 1 round, S-Box limit at least increase 1
   while True:
      cnfdoc = FormulateCNF(round, SBoxLimit, MatsuiBounds)
      Runtime, Solution = CallSolver(cnfdoc, solver_choice)
      if (Runtime == -1):
         SBoxLimit += 1
      else:
         # Stop Timer
         time_end = time.time()
         # Record round information
         DiffActiveSbox[round] = SBoxLimit   # Update active S-box limit
         SolutionSet.append([round, SBoxLimit, Runtime, Solution])   # Record solution
         with open("./runlog_" + solver_choice + '_' + type_matsui + '_' + str(target_round), "a") as log:
            log.write("Round" + "{0:<3}".format(str(round)+':') + " Minimum Active SBox = " + "{0:<2}".format(str(SBoxLimit)) + "; Solving Cost = " + "{0:1.6f}".format(Runtime) + "; Cumulative Total Cost = " + "{0:1.6f}".format(time_end-total_start) + "\n")
         break
