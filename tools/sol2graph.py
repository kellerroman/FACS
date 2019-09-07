#!/usr/bin/env python3

from pyx import *

file = open("sol.dat","r")

frist = file.readline()

first = frist.split()

nCells = int(first[0])
nPnts  = int(first[1])
Q_DIM  = int(first[2])

print(nCells,nPnts,Q_DIM)
pnts = []
cells = []
neigh = []
for i in range(nPnts):
    line = (file.readline())
    line = line.split(":")[1]
    line = line.split()
    #print(line)
    pnts.append([float(line[0])*20.,float(line[1])*20.])

for i in range(nCells):
    line = (file.readline())
    line = line.split(":")[1]
    line = line.split()
    #print(line)
    cells.append([int(line[0])-1,int(line[1])-1,int(line[2])-1,int(line[3])-1])

for i in range(nCells):
    line = (file.readline())
    line = line.split(":")[1]
    line = line.split()
    #print(line)
    neigh.append([int(line[0])-1,int(line[1])-1,int(line[2])-1,int(line[3])-1])
file.close()

#print(pnts)
#print(cells)

c = canvas.canvas()
for i in range(nCells):
    for j in range(4):
        p1 = cells[i][j]
        p2 = cells[i][(j+1) % 4]
#        print(p1," => ",p2)
        c.stroke(path.line(pnts[p1][0],pnts[p1][1], pnts[p2][0], pnts[p2][1]),[style.linewidth(0.0001)])


    tmp = pnts[cells[i][0]][1] + pnts[cells[i][3]][1]
    tmp *= 0.5
    t = str(i+1)
    t+="("+str(neigh[i][0]+1)
    t+=","+str(neigh[i][1]+1)
    t+=","+str(neigh[i][2]+1)
    t+=","+str(neigh[i][3]+1)+")"
    c.text(pnts[cells[i][0]][0],tmp, t, [text.halign.boxleft])
    # c.text(pnts[cells[i][0]][0],tmp, r"text", [text.halign.boxleft])

#c.stroke(path.line(0, 0, 1, 0))
#c.stroke(path.line(1, 0, 1, 1))

c.writeEPSfile("path")
#c.writeSVGfile("path")
