import optparse
import os

def main():
    parser = optparse.OptionParser()
    parser.add_option("-p", "--psicov", dest="psc", help="path to psicov alignment file", metavar="STRING")
    parser.add_option("-c", "--combs", dest="combs", help="path to combs out file", metavar="STRING" )

    (options, args) = parser.parse_args()

    psc = options.psc
    combs_out = options.combs

    with open(psc, 'r') as pscfile:
        combs_seq = pscfile.readline()

    name = os.path.basename(psc).split(".")[0]

    with open(combs_out,"w") as o:
        o.write(">%s %i\n%s\n"%(name, len(combs_seq),combs_seq))


    print("Successfully extracted combs sequence from psicov file for {0}".format(name))



if __name__ == '__main__':
    main()