
from re import sub

inp = open('patterns-c.txt', 'r')

data = []

for (i, line) in enumerate(inp):
    line = line.strip()
    data.append(line.split(' '))

n = len(data)

print('<table>')

ncol = 7

for i in range((n+ncol-1) // ncol):
    print('<tr>')
    for j in range(ncol):
        if (ncol*i+j < n):
            print('<td>%s</td>' % data[ncol*i+j][0])
    print('</tr>')
    print('<tr>')
    for j in range(ncol):
        if (ncol*i+j < n):
            print('<td>%s</td>' % data[ncol*i+j][1])
    print('<tr>')
    for j in range(ncol):
        if (ncol*i+j < n):
            print('<td><img src="img/%s.png"/></td>' % sub(',', '', data[ncol*i+j][0]))
    print('</tr>')

print('</table>')
