import sys, os
from operator import itemgetter


def get_pdb(pdb):

        pdbf = '/salilab/park2/database/pdb/divided/%s/pdb%s.ent.gz' % (pdb[1:3].lower(), pdb.lower())

        os.system("cp %s ." % (pdbf,))
        # untar
        os.system("gunzip pdb*.ent.gz")
        os.system('mv pdb%s.ent %s.pdb' % (pdb.lower(),pdb))


def get_chains(pdb, chain, ligand=0):
        data = open(pdb+'.pdb')
        D = data.readlines()
        data.close()

	AA = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
	      'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
   	      'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    	      'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
	Seq = []
        for d in D:
                if 'ATOM'==d[:4]:
                        if d[21]==chain and d[17:20]!= 'HOH':
                                res, rsid = d[17:20],int(d[22:26])
				try: Seq.append((rsid,AA[res]))
				except KeyError: pass
	return ''.join([i[1] for i in sorted(list(set(Seq)), key=itemgetter(0))])


protein= sys.argv[1]
hhfile = sys.argv[2]
pdbids = sys.argv[3:]

data = open(hhfile)
D = data.read()
data.close()

print
print
for pdb in pdbids:
	
	query = ''
	templ = ''
	chain = ''
	qstart,qend = 10000,0
	tstart,tend = 10000,0
	for d in D.split('>')[1:]:
		if pdb==d[:4]:
			chain = d[5]
			d = d.split('\n')[5:-4]
			for l in xrange(0,len(d),10):
				query += d[l].split()[3]
				templ += d[l+4].split()[3]

				if int(d[l].split()[2])<qstart: qstart=int(d[l].split()[2])
				if int(d[l].split()[4])>qend: qend=int(d[l].split()[4])

				if int(d[l+4].split()[2])<tstart: tstart=int(d[l+4].split()[2])
				if int(d[l+4].split()[4])>tend: tend=int(d[l+4].split()[4])

		
	get_pdb(pdb)
	#print get_chains(pdb,chain)	


	print '>P1;%s_%i_%i' % (protein,qstart,qend)
	print 'sequence:%s_%i_%i:  : :  : ::::' % (protein,qstart,qend)
	print '\n'.join([query[x:x+100] for x in xrange(0,len(query),100)])+'*'
	print
	print '>P1;%s' % (pdb,)
	print 'structure:%s:  %i: %s: %i :%s ::::' % (pdb,tstart,chain,tend,chain)
	print '\n'.join([templ[x:x+100] for x in xrange(0,len(templ),100)])+'*'
	print
	print
	

