#!/usr/bin/env python3

from pyx import *

file = open("sol.dat","r")

line = (file.readline())
frist = file.readline()

first = frist.split()

nCells = int(first[0])
nPnts  = int(first[1])
Q_DIM  = int(first[2])

print(nCells,nPnts,Q_DIM)
pnts = []
center = []
cells = []
neigh = []
line = (file.readline())
for i in range(nPnts):
    line = (file.readline())
    line = line.split(":")[1]
    line = line.split()
    #print(line)
    pnts.append([float(line[0])*20.,float(line[1])*20.])

line = (file.readline())
for i in range(nCells):
    line = (file.readline())
    line = line.split(":")[1]
    line = line.split()
    #print(line)
    cells.append([int(line[0])-1,int(line[1])-1,int(line[2])-1,int(line[3])-1])

line = (file.readline())
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
    cx = 0
    cy = 0
    for j in range(4):
        p1 = cells[i][j]
        p2 = cells[i][(j+1) % 4]
#        print(p1," => ",p2)
        c.stroke(path.line(pnts[p1][0],pnts[p1][1], pnts[p2][0], pnts[p2][1]),[style.linewidth(0.001)])
        cx += pnts[p1][0]
        cy += pnts[p1][1]

    cx *= 0.25
    cy *= 0.25
    center.append([cx,cy])

for i in range(nCells):
    for j in range(4):
        n = neigh[i][j]
        if n > -1:
            nx = center[n][0]
            ny = center[n][1]
            mx = center[i][0]
            my = center[i][1]
            d = ((nx-mx)**2+(ny-my)**2)**(0.5) * 0.05
            if mx > nx:
                mx -= d
                nx += d
            elif nx > mx:
                mx += d
                nx -= d
            if my > ny:
                my -= d
                ny += d
            elif ny > my:
                my += d
                ny -= d
            c.stroke(path.line(mx,my,nx,ny)
                ,[style.linewidth(0.0001), color.rgb.red, deco.earrow])

for i in range(nCells):
    t = str(i+1)
    # t+="("+str(neigh[i][0]+1)
    # t+=","+str(neigh[i][1]+1)
    # t+=","+str(neigh[i][2]+1)
    # t+=","+str(neigh[i][3]+1)+")"
    c.text(center[i][0], center[i][1], t, [text.halign.boxcenter, text.halign.flushcenter, text.valign.middle])

c.writeEPSfile("path")
#c.writeSVGfile("path")
