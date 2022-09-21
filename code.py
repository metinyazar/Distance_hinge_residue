from prody import *

## pdb_id : 4 letter code of pdb files
## res_pos : Postion of residue
## chain_id: chain of the residue

def GNM_calc_distance(pdb_id,res_pos,chain_id):
    pdb = parsePDB(pdb_id + '.pdb')
    calphas = pdb.select('calpha')
    gnm = GNM()
    gnm.buildKirchhoff(calphas)
    gnm.calcModes()
    resnums = calphas.getResnums()
    chains = calphas.getChids()
    mode_hinges = calcHinges(gnm[0])
    hinges_res = resnums[mode_hinges]
    hinges_chain = chains[mode_hinges]
    atom1 = pdb.select('calpha and chain ' + chain_id + ' and resnum ' + res_pos)
    data = []
    for i,j in zip(hinges_res,hinges_chain):
        atom2 = pdb.select('calpha and chain '+str(j)+ ' and resnum '+str(abs(i)))
        data.append(list(calcDistance(atom1,atom2)))
    data2 = [j[0] for j in data]
    return min(data2)
