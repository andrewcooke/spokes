
from re import sub

inp = open('patterns.txt', 'r')

data = []

for (i, line) in enumerate(inp):
    line = line.strip()
    data.append(line.split(' '))

n = len(data)

print('<table>')

ncol = 4

for i in range((n+1) // ncol):
    print('<tr>')
    for j in range(ncol):
        if i == j == 0:
            print('<td>Name</td>')
            print('<td>Length</td>')
        else:
            print('<td>%s</td>' % data[ncol*i+j-1][0])
            print('<td>%s</td>' % data[ncol*i+j-1][1])
    print('</tr>')
    print('<tr>')
    for j in range(ncol):
        if i == j == 0:
            print('<td colspan="2">Image of the pattern (highlighted in red) used in a typical wheel (20, 32 or 36 spokes).</td>')
        else:
            print('<td colspan="2"><img src="img/%s.png"/></td>' % sub(',', '', data[ncol*i+j-1][0]))
    print('</tr>')
    print('<tr>')
    for j in range(ncol):
        if i == j == 0:
            print('<td colspan="2">Common names</td>')
        else:
            print('<td colspan="2"></td>')
    print('</tr>')


print('</table>')
