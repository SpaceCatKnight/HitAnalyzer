import ROOT as rt
 
# open tree file
file = rt.TFile("TT1.root",'READ')
tree = file.Get("demo/tree")

N = tree.GetEntries()
for i in xrange(N):
    print "Working on event " ,i
    tree.GetEntry(i)
    print "This event has %i AK4 jets" %tree.nJets
    for j in range(0,tree.nJets):
       print "For jet " ,j
       print "Jet pt = " , tree.jet_pt[j]
       print "Jet eta = " , tree.jet_eta[j]
       print "Jet phi = " , tree.jet_phi[j]
       print "Jet mass = " , tree.jet_mass[j]; print "";
    

 