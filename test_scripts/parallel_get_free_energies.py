
import multiprocessing as mp
import glob

PATH = "/home/jeanellg/biodata/maize/rawdata/GENOME/contigs/"

def get_lowest(path):
        data = glob.glob(path + '*.ct')
        structures = dict()
        for i in data:
                data_file = open(i, 'r')
                energies = data_file.readline().strip().split()
                fin_i = energies.index('=') + 1
                if "[initially" in energies:
                        init_i = energies.index('[initially') + 1
                        init = float(energies[init_i].replace(']',''))
                else:
                        init = "not in file"
                fin = float(energies[fin_i])
                name = i.split('/')[-1]
                structures[name] = (init,fin)
                data_file.close()
        lowest = structures.items()
        lowest.sort(key = lambda x:(abs(x[1][1]),x[0]))
        ret = []
        for x in lowest:
                if (str(x[1][1]) == str(lowest[0][1][1])):
                        ret.append(x)
                else:
                        break
        return ret

if __name__ == "__main__":
        paths = glob.glob(PATH + "*/")
        if (len(paths) < 10):
                cores = len(paths)
        else:
                cores = 10
        pool = mp.Pool(cores)
        results = pool.map(get_lowest,paths)
        w_file = open('all_fin_free_energies.txt','w')
        w_file.write("file\tinit_energy\tfin_energy\n")
        print("writing to file")
        for data in results:
                for x in data:
                        w_file.write('{}\t{}\t{}\n'.format(x[0],x[1][0],x[1][1]))
        w_file.close()
