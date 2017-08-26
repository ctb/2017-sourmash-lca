import json

info = {}
info['version'] = 1

db = []
info['dblist'] = db

k31 = {}
k31['ksize'] = 31
k31['lca_db'] = 'genbank-k31.lca'
k31['scaled'] = 10000
k31['nodes'] = 'genbank/nodes.dmp'
k31['names'] = 'genbank/names.dmp'
db.append(k31)

k51 = {}
k51['ksize'] = 51
k51['lca_db'] = 'genbank-k51.lca'
k51['scaled'] = 10000
k51['nodes'] = 'genbank/nodes.dmp'
k51['names'] = 'genbank/names.dmp'
db.append(k51)


k21 = {}
k21['ksize'] = 21
k21['lca_db'] = 'genbank-k21.lca'
k21['scaled'] = 10000
k21['nodes'] = 'genbank/nodes.dmp'
k21['names'] = 'genbank/names.dmp'
db.append(k21)

json.dump(info, open('genbank.lca.json', 'wt'))
