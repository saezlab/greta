import mudata, sys

m = mudata.read(sys.argv[1])
m.mod['scRNA'] = m.mod['rna']
del m.mod['rna']
m.mod['scATAC'] = m.mod['atac']
del m.mod['atac']
m.mod['scATAC'].var_names = m.mod['scATAC'].var_names.str.replace('-', ':', 1)
m.write(sys.argv[2])
