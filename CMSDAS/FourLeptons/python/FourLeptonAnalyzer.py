from CMSDAS.FourLeptons.Analyzer import *

class FourLeptonAnalyzer(Analyzer):
    def __init__(self):
        super(FourLeptonAnalyzer,self).__init__()






    #####CHANGE THIS METHOD TO CHANGE MUON ID######    
    def muonID(self,muon,vertex):
        if muon.pt()<5 or abs(muon.eta())>2.4:
            return False
        if muon.innerTrack().dxy( vertex.position() )>0.02:
            return False
        if muon.innerTrack().dz( vertex.position() )>0.2:
            return False
        if not (muon.isPFMuon() and \
               ( muon.isGlobalMuon() or muon.isTrackerMuon() )):
            return False
        # muon ISO variable
        if (muon.chargedHadronIso()+max(0.0,muon.photonIso()+muon.neutralHadronIso()-0.5*muon.puChargedHadronIso()))/muon.pt()>0.6:
            return False
        # muon SIP variable
        if (muon.dB(2)/muon.edB(2))>8:
            return False

        return True
    

    
    #####CHANGE THIS METHOD TO CHANGE ELECTRON ID######    
    def electronID(self,electron,vertex):
        if electron.pt()<7 or abs(electron.eta())>2.5:
            return False


        mvaRegions = [{'ptMin':0,'ptMax':10, 'etaMin':0.0, 'etaMax':0.8,'mva':0.47},\
                      {'ptMin':0,'ptMax':10, 'etaMin':0.8 ,'etaMax':1.479,'mva':0.004},\
                      {'ptMin':0,'ptMax':10, 'etaMin':1.479, 'etaMax':3.0,'mva':0.295},\
                      {'ptMin':10,'ptMax':99999999, 'etaMin':0.0, 'etaMax':0.8,'mva':-0.34},\
                      {'ptMin':10,'ptMax':99999999, 'etaMin':0.8, 'etaMax':1.479,'mva':-0.65},\
                      {'ptMin':10,'ptMax':99999999, 'etaMin':1.479, 'etaMax':3.0,'mva':0.6}]
        ID=False 
        for element in mvaRegions:
            if electron.pt()>= element['ptMin'] and \
               electron.pt()< element['ptMax'] and \
               abs(electron.superCluster().eta())>=element['etaMin'] and \
               abs(electron.superCluster().eta())<element['etaMax'] and \
               electron.electronID("mvaNonTrigV0")> element['mva']: 
                ID=True
        if not ID:
            return False

        if electron.gsfTrack().trackerExpectedHitsInner().numberOfHits()>1:
            return False

        # electron ISO variable
        if (electron.chargedHadronIso()+max(0.0,electron.photonIso()+electron.neutralHadronIso()-0.5*electron.puChargedHadronIso()))/electron.pt()>0.6:
            return False
        # electron SIP variable
        if (electron.dB(2)/electron.edB(2))>8:
            return False


        return True


    def analyze(self,box):

        #####START FROM A bOX CONTAINING SELECTED MUONS AND ELECTRONS and MAKE
        #FOUR LEPTON CANDIDATES
        
        
        #Now check if there are at least four leptons:
        box.leptons=set(box.selectedMuons+box.selectedElectrons)
        if len(box.leptons)<4:
            return False

        #Now create Z candidates and apply cuts:
        box.zcandidates=[]

        for l1,l2 in itertools.combinations(box.leptons,2):
            #ask that lepton ID passed for the Z1 ONLY
            #for data driven estimation
            #they need to have same flavour and OS
            if abs(l1.pdgId()) != abs(l2.pdgId()):
                continue
            if l1.charge() +l2.charge() !=0:
                continue
            #now create a di lepton object and check mass
            z=DiObject(l1,l2)
            if not (z.mass()>12 and z.mass()<120):
                continue
            box.zcandidates.append(DiObject(l1,l2))
        # OK if there are more than one Z candidates
        #pick the one with the best mass.
        if len(box.zcandidates)==0:
            return False



        sortedZs=sorted(box.zcandidates,key=lambda x: abs(x.mass()-91.118))
        box.Z1 = sortedZs[0]
        
        #now remove the used leptons from the list and make Z2 pairs
        box.leptons.remove(box.Z1.l1)        
        box.leptons.remove(box.Z1.l2)


        #now the same thing with the second
        box.zcandidates2=[]

        for l1,l2 in itertools.combinations(box.leptons,2):
            #they need to have same flavour and OS
            if abs(l1.pdgId()) != abs(l2.pdgId()):
                continue
            if l1.charge() +l2.charge() !=0:
                continue
            #now create a di lepton object and check mass
            z=DiObject(l1,l2)
            if not (z.mass()>4 and z.mass()<120):
                continue
            box.zcandidates2.append(DiObject(l1,l2))
        # OK if there are more than one Z candidates
       #pick the one with the highest lepton pt sum
        if len(box.zcandidates2)==0:
            return False


        
        sortedZ2s=sorted(box.zcandidates2,key=lambda x: x.l1.pt()+x.l2.pt(),reverse=True)
        box.Z2 = sortedZ2s[0]

        #kill the candidate if a OS pair has mll<4 GeV
        for l1,l2 in itertools.combinations([box.Z1.l1,box.Z1.l2,box.Z2.l1,box.Z2.l2],2):
            ll =DiObject(l1,l2)
            if (l1.charge()+l2.charge()) ==0:
                if ll.mass()<4 :
                    return False


        
        #create the ZZ
        box.ZZ = DiObject(box.Z1,box.Z2)

        # compute 5 angles [costhetastar, costheta1, costheta2, Phi, Phi1]
        angles = computeAngles(box.ZZ.l1.l1, box.ZZ.l1.l2, box.ZZ.l2.l1, box.ZZ.l2.l2)
        box.costhetastar = angles[0]
        box.costheta1 = angles[1]
        box.costheta2 = angles[2]
        
        return True



    def declareHistos(self):
        super(FourLeptonAnalyzer,self).declareHistos()
        ###ADD YOUR HISTOGRAMS AFTER THIS LINE AS AbOVE#####
        self.declareHisto('mass',30,70,150,"m_{4l} [GeV]")
        self.declareHisto('massFull',100,70,570,"m_{4l} [GeV]")
        self.declareHisto('massZ1',20,12,120,"m_{Z1} [GeV]")
        self.declareHisto('massZ2',20,4,74,"m_{Z2} [GeV]")

        # anlges
        self.declareHisto('costhetastar',20,-1,+1,"cos(#theta^{*})")
        self.declareHisto('costheta1',20,-1,+1,"cos(#theta_{1})")
        self.declareHisto('costheta2',20,-1,+1,"cos(#theta_{2})")


    def fillHistos(self,box,sample,weight = 1):
        super(FourLeptonAnalyzer,self).fillHistos(box,sample,weight)

        self.fillHisto('mass',sample,box.ZZ.mass(),weight)        
        self.fillHisto('massFull',sample,box.ZZ.mass(),weight)        
        self.fillHisto('massZ1',sample,box.ZZ.l1.mass(),weight)        
        self.fillHisto('massZ2',sample,box.ZZ.l2.mass(),weight)        

        # anlges
        self.fillHisto('costhetastar',sample,box.costhetastar,weight)
        self.fillHisto('costheta1',sample,box.costheta1,weight)
        self.fillHisto('costheta2',sample,box.costheta2,weight)



# compute 5 angles for given 4 leptons, according to http://arxiv.org/pdf/1208.4018v2.pdf
def computeAngles(l11, l12, l21, l22):
                
    # get p4 for Z1 and Z2
    p4M11 = ROOT.TLorentzVector(l11.px(), l11.py(), l11.pz(), l11.energy())
    p4M12 = ROOT.TLorentzVector(l12.px(), l12.py(), l12.pz(), l12.energy())
    p4M21 = ROOT.TLorentzVector(l21.px(), l21.py(), l21.pz(), l21.energy())
    p4M22 = ROOT.TLorentzVector(l22.px(), l22.py(), l22.pz(), l22.energy())

    # for OS pairs: lep1 must be the negative one; for SS pairs: use random deterministic convention
    if ((l11.pdgId()*l12.pdgId() < 0 and l11.pdgId()<0) or (l11.pdgId()*l12.pdgId() > 0 and l11.phi()<l12.phi())):
	p4M11, p4M12 = p4M12, p4M11
    if ((l21.pdgId()*l22.pdgId() < 0 and l21.pdgId()<0) or (l21.pdgId()*l22.pdgId() > 0 and l21.phi()<l22.phi())):
        p4M21, p4M22 = p4M22, p4M21

    # get p4 for Z1 and Z2
    p4Z1 = p4M11 + p4M12
    p4Z2 = p4M21 + p4M22

    ### computation of cos(thetaStar) ###
    # prepare p4 for H, use it to boost p4Z1 and p4Z2
    p4H = p4Z1 + p4Z2
    boostX = -(p4H.BoostVector())
    thep4Z1inXFrame = p4Z1
    thep4Z2inXFrame = p4Z2
    thep4Z1inXFrame.Boost( boostX )
    thep4Z2inXFrame.Boost( boostX )

    # compute cos(thetaStar) from p3 vectors of Z1 and Z2 in X ref. frame
    theZ1X_p3 = ROOT.TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() )
    theZ2X_p3 = ROOT.TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() )
    costhetastar = theZ1X_p3.CosTheta()

    ### computation of cos(theta1) ###
    # boost lep. p4 to the Z1 ref. frame
    boostV1 = -(p4Z1.BoostVector())
    p4M11_BV1 = p4M11
    p4M21_BV1 = p4M21
    p4M22_BV1 = p4M22
    p4M11_BV1.Boost( boostV1 )
    p4M21_BV1.Boost( boostV1 )
    p4M22_BV1.Boost( boostV1 )

    # compute cos(theta1) from q11 and q2 in Z1 ref. frame 
    q11 = p4M11_BV1.Vect()
    q2 = (p4M21_BV1 + p4M22_BV1).Vect()
    costheta1 = (-1)*( q2.Dot( q11 ) )/( q2.Mag() * q11.Mag() );

    ### computation of cos(theta2) ###
    # boost lep. p4 to the Z2 ref. frame
    boostV2 = -(p4Z2.BoostVector());
    p4M11_BV2 = p4M11
    p4M21_BV2 = p4M21
    p4M12_BV2 = p4M12
    p4M11_BV2.Boost( boostV2 )
    p4M21_BV2.Boost( boostV2 )
    p4M12_BV2.Boost( boostV2 )

    # compute cos(theta2) from q21 and q1 in Z2 ref. frame
    q21 = p4M21_BV2.Vect()
    q1 = (p4M11_BV2 + p4M12_BV2).Vect()
    costheta2 = (-1)*( q1.Dot( q21 ) )/( q1.Mag() * q21.Mag() );

    # return list of angles [costhetastar, costheta1, costheta2, Phi, Phi1]
    angles = [costhetastar, costheta1, costheta2, 0, 0]
    return angles


        
        
