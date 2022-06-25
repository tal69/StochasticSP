"""  Create a LaTex graph that demonstrate the accuracy of piece-wise approximation """

import Single_SP
import numpy as np
C = 30
p = 0.5
small_ticks = 91
jumps = 9
print("starting...")
Lambda = np.array([C*p*i/60 for i in range(small_ticks)])
R = np.array([Single_SP.markov_chain_sp(C,Lambda[i],p)[0] for i in range(small_ticks)])

Lambda11 = np.array([Lambda[i*jumps] for i in range(11)])
R11 = np.array([R[i*jumps] for i in range(11)])

print("""
\\begin{tikzpicture}
\\begin{axis}[
    xlabel={$\lambda$},
    ylabel={Number of rejections},""")
print(f"    xmin={Lambda[0]}, xmax={Lambda[-1]},")
print(f"    ymin={R[0]}, ymax={np.ceil(R[-1])},")

    # xtick={35,40,45,50},
    # ytick={3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900},
    # legend pos=south east,
    # ymajorgrids=true,
    # grid style=dashed,
print("]")


print("""
\\addplot[
    color=red,
    mark=circle,
    ]
    coordinates {""")
for i in range(small_ticks):
    print(f"({Lambda[i]:.5f},{R[i]:.5f})", end="")
print("};")

print("""
\\addplot[
    color=blue,
    mark=circle,
    ]
    coordinates {""")
for i in range(11):
    print(f"({Lambda11[i]:.5f},{R11[i]:.5f})", end="")
print("};")

print("""
\\end{axis}
\\end{tikzpicture}""")


"""  Find maximum error """

print(max([np.interp(Lambda[i], Lambda11,  R11) - R[i]  for i in range(91)]))
print(Lambda[np.argmax([np.interp(Lambda[i], Lambda11,  R11) - R[i]  for i in range(91)])])
