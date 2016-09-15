
from re import sub

inp = open('patterns.txt', 'r')

data = []

for (i, line) in enumerate(inp):
    line = line.strip()
    data.append(line.split(' '))

n = len(data)

print('<table>')

ncol = 4

for i in range(n // ncol):
    print('<tr>')
    for j in range(ncol):
        print('<td>%s</td>' % data[ncol*i+j][0])
        print('<td>%s</td>' % data[ncol*i+j][1])
    print('</tr>')
    print('<tr>')
    for j in range(ncol):
        print('<td colspan="2"><img src="img/%s.png"/></td>' % sub(',', '', data[ncol*i+j][0]))
    print('</tr>')
    print('<tr>')
    for j in range(ncol):
        print('<td colspan="2"></td>')
    print('</tr>')


print('</table>')
