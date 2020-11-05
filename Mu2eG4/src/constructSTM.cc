//
// Free function to create the Stopping Target Monitor (STM)
//
// Author: Anthony Palladino
//
// Notes:
//
// The initial implementaion is described in Mu2e Document XXXX


// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e includes.
#include "Mu2eG4/inc/constructSTM.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/G4GeometryOptions.hh"
#include "STMGeom/inc/STM.hh"
#include "STMGeom/inc/PermanentMagnet.hh"
#include "STMGeom/inc/SupportTable.hh"
#include "STMGeom/inc/TransportPipe.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "GeometryService/inc/VirtualDetector.hh"
#include "DataProducts/inc/VirtualDetectorId.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/SensitiveDetectorHelper.hh"
#include "GeomPrimitives/inc/PolyconsParams.hh"

// G4 includes
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Torus.hh"
#include "G4Para.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"

using namespace std;

namespace mu2e {

  void constructSTM(const SimpleConfig& _config){

    STM const & stmgh = *(GeomHandle<STM>());

    PermanentMagnet const & pSTMMagnetParams               = *stmgh.getSTMMagnetPtr();
    TransportPipe   const & pSTMTransportPipeParams        = *stmgh.getSTMTransportPipePtr();
    STMCollimator   const & pSTMFOVCollimatorParams        = *stmgh.getSTMFOVCollimatorPtr();
    SupportTable    const & pSTMMagnetSupportTableParams   = *stmgh.getSTMMagnetSupportTablePtr();
    STMCollimator   const & pSTMSSCollimatorParams         = *stmgh.getSTMSSCollimatorPtr();
	STMShielding   const & pSTMShieldingParams         = *stmgh.getSTMShieldingPtr();
    SupportTable    const & pSTMDetectorSupportTableParams = *stmgh.getSTMDetectorSupportTablePtr();
	STMAbsorber     const & pSTMAbsorberParams             = *stmgh.getSTMAbsorberPtr();//wyq
	STMBasorber     const & pSTMBasorberParams             = *stmgh.getSTMBasorberPtr();//wyq
//    GeDetector      const & pSTMDetector1Params            = *stmgh.getSTMDetector1Ptr();
    GeDetector      const & pSTMDetector2Params            = *stmgh.getSTMDetector2Ptr();
    ShieldPipe      const & pSTMShieldPipeParams           = *stmgh.getSTMShieldPipePtr();

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "stm", "stm");

    const bool STMisVisible        = geomOptions->isVisible("stm");
    const bool STMisSolid          = geomOptions->isSolid("stm");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("stm");
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("stm");
    const bool placePV             = geomOptions->placePV("stm");
    int  const verbosityLevel      = _config.getInt("stm.verbosityLevel", 0);
    const double stmZAllowed=_config.getDouble("stm.z.allowed");
    const G4ThreeVector zeroVector(0.,0.,0.);
		const double stmSSCollxshift=_config.getDouble("stm.SScollimator.xshift");
		const double stmSSCollyshift=_config.getDouble("stm.SScollimator.yshift");
		const int stmSSCollShape=_config.getInt("stm.SScollimator.shape",1);
G4ThreeVector stmSSCollShift(stmSSCollxshift,stmSSCollyshift,0.0);
//Read in parameters from configuration
const double STMShield_Ttop_Poly=_config.getDouble("STMShield_Ttop_Poly");
const double STMShield_Ttop_Pb1=_config.getDouble("STMShield_Ttop_Pb1");
const double STMShield_Ttop_Pb2=_config.getDouble("STMShield_Ttop_Pb2");
const double STMShield_Ttop_Cu=_config.getDouble("STMShield_Ttop_Cu");
const double STMShield_Ttop_Al=_config.getDouble("STMShield_Ttop_Al");
const double STMShield_Tbot_Al=_config.getDouble("STMShield_Tbot_Al");
const double STMShield_Ttop_Ycrack=_config.getDouble("STMShield_Ttop_Ycrack",0.0);
const double STMShield_X1=_config.getDouble("STMShield_X1");
const double STMShield_X2=_config.getDouble("STMShield_X2");
const double STMShield_Y1=_config.getDouble("STMShield_Y1");
const double STMShield_Zgap=_config.getDouble("STMShield_Zgap");
const double STMDet2HPGe_Can_R=_config.getDouble("STMDet2HPGe_Can_R");
const double STMDet2HPGe_Can_thick=_config.getDouble("STMDet2HPGe_Can_thick");
const double STMDet2HPGe_Can_halfL=_config.getDouble("STMDet2HPGe_Can_halfL");
const double STMDet2HPGe_Can_gap_top=_config.getDouble("STMDet2HPGe_Can_gap_top");
const double STMDet2HPGe_Cup0_gap_bar=_config.getDouble("STMDet2HPGe_Cup0_gap_bar");
const double STMDet2HPGe_Cup0_gap_bot=_config.getDouble("STMDet2HPGe_Cup0_gap_bot");
const double STMDet2HPGe_Cup0_thick_bar=_config.getDouble("STMDet2HPGe_Cup0_thick_bar");
const double STMDet2HPGe_Cup0_thick_top=_config.getDouble("STMDet2HPGe_Cup0_thick_top");
const double STMDet2HPGe_Cup0_thick_bot=_config.getDouble("STMDet2HPGe_Cup0_thick_bot");
const double STMDet2HPGe_R_out_corner=_config.getDouble("STMDet2HPGe_R_out_corner");
const double STMDet2HPGe_R_out=_config.getDouble("STMDet2HPGe_R_out");
const double STMDet2HPGe_L_out=_config.getDouble("STMDet2HPGe_L_out");
const double STMDet2HPGe_R_in_corner=_config.getDouble("STMDet2HPGe_R_in_corner");
const double STMDet2HPGe_R_in=_config.getDouble("STMDet2HPGe_R_in");
const double STMDet2HPGe_L_in=_config.getDouble("STMDet2HPGe_L_in");
const double STMDet2HPGe_xOffset=_config.getDouble("STMDet2HPGe_xOffset");
//>LaBr3
const double STMShield_LaBr3_halfWidth=_config.getDouble("STMShield_LaBr3_halfWidth");
const double STMShield_LaBr3_halfT_Al=_config.getDouble("STMShield_LaBr3_halfT_Al");
const double STMShield_LaBr3_gap_Al=_config.getDouble("STMShield_LaBr3_gap_Al");
const double STMShield_LaBr3_gap_Pb=_config.getDouble("STMShield_LaBr3_gap_Pb");
const double STMShield_LaBr3_halfT_Cu=_config.getDouble("STMShield_LaBr3_halfT_Cu");
const double STMShield_LaBr3_halfT_Pb=_config.getDouble("STMShield_LaBr3_halfT_Pb");
const double STMShield_LaBr3_Hole_R=_config.getDouble("STMShield_LaBr3_Hole_R");
//Shield btw. HPGe and LaBr3
const double STMShield_halfLengthI=_config.getDouble("STMShield_halfLengthI");
const double STMShield_halfTI_Pb=_config.getDouble("STMShield_halfTI_Pb");
const double STMShield_halfTI_Al=_config.getDouble("STMShield_halfTI_Al");
const double STMShield_halfTI_Cu=_config.getDouble("STMShield_halfTI_Cu");
//Parameters for Det1(LaBr3) and Det2(HPGe)
const double STMDet1CupLaBr3_R=_config.getDouble("STMDet1CupLaBr3_R");
const double STMDet1CupLaBr3_thick=_config.getDouble("STMDet1CupLaBr3_thick");
const double STMDet1CupLaBr3_halfLength=_config.getDouble("STMDet1CupLaBr3_halfLength");
const double STMDet1CupLaBr3_Upstream=_config.getDouble("STMDet1CupLaBr3_Upstream");
const double STMDet1LaBr3_R=_config.getDouble("STMDet1LaBr3_R");
const double STMDet1LaBr3_halfLength=_config.getDouble("STMDet1LaBr3_halfLength");
const double STMDet1LaBr3_gap_Upstream=_config.getDouble("STMDet1LaBr3_gap_Upstream");
const double STMDet1LaBr3_xOffset=_config.getDouble("STMDet1LaBr3_xOffset");

double STMShield_Z1=stmZAllowed-pSTMSSCollimatorParams.halfLength()*2.0;//901.7;
double STMShield_Z2=STMShield_Z1-STMShield_Zgap-(STMShield_Tbot_Al+2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Cu);
double STMShield_Theta=-atan((STMShield_X2-STMShield_X1)/2.0/STMShield_Z1);//-0.41757;
double STMShield_right_Theta=atan((STMShield_X2-STMShield_X1)/STMShield_Z1);//top
double STMShield_thick_front=2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Cu+STMShield_Tbot_Al;
double STMShield_thick_bot=2*STMShield_Ttop_Poly+2*STMShield_Ttop_Pb2+STMShield_Ttop_Cu+STMShield_Tbot_Al;
double STMShield_Ttop_CenterX=(STMShield_X2+STMShield_X1)/4.0-STMShield_X1/2.0;//dis from center to straight side
    if ( verbosityLevel > 0) {
       cout << __func__ << " STM verbosityLevel    : " << verbosityLevel  << endl;
    }

    // Access to the G4HelperService.
    G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));

    //For now just make the mother HallAir because the STM transpport pipe goes
    //into the End Cap Shielding mother
    //TODO: create STM Mother volume and make that the parent instead of HallAir
    VolumeInfo const & parentInfo = _helper->locateVolInfo("HallAir");
    G4ThreeVector parentCenterInMu2e = parentInfo.centerInMu2e();

    // Fetch DS geom. object
    GeomHandle<DetectorSolenoid> ds;
    const CLHEP::Hep3Vector &dsP( ds->position() );

    GeomHandle<CosmicRayShield> CRS;
    std::vector<double> crvd_halflengths = CRS->getSectorHalfLengths("D");
    CLHEP::Hep3Vector   crvd_position    = CRS->getSectorPosition("D");
    const double z_crv_max = crvd_position.z() + crvd_halflengths[2];

    // Create a reference position (most things in the STM geometry will be defined w.r.t. this position)
    // Our reference z is the downstream edge of the CRV
    const G4ThreeVector stmReferencePositionInMu2e(dsP.x(),
                                                    0.0,
                                                    z_crv_max );
    const G4ThreeVector stmReferencePositionInParent = stmReferencePositionInMu2e - parentCenterInMu2e;


    //===================== Sweeper Magnet ==========================

    G4ThreeVector stmMagnetPositionInMu2e   = pSTMMagnetParams.originInMu2e();
    G4ThreeVector stmMagnetPositionInParent = pSTMMagnetParams.originInMu2e() - parentCenterInMu2e;

    // Make the magnet
    G4Box* boxMagnet     = new G4Box("boxMagnet",     pSTMMagnetParams.xHalfLength(),     pSTMMagnetParams.yHalfLength(),     pSTMMagnetParams.zHalfLength());

    // Make the rectangular window (make the box that gets subtracted just a bit longer to be sure there are no edge effects)
    const double stmMagnetHoleHalfLengths[3] = {pSTMMagnetParams.xHoleHalfLength(),
                                                pSTMMagnetParams.yHoleHalfLength(),
                                                pSTMMagnetParams.zHalfLength()     };
    G4Box* boxMagnetHole = new G4Box("boxMagnetHole", stmMagnetHoleHalfLengths[0]+pSTMShieldPipeParams.linerWidth(), stmMagnetHoleHalfLengths[1]+pSTMShieldPipeParams.linerWidth(), stmMagnetHoleHalfLengths[2]+1.0);
    VolumeInfo stmMagnet;
    stmMagnet.name = "stmMagnet";
    stmMagnet.solid = new G4SubtractionSolid(stmMagnet.name,boxMagnet,boxMagnetHole,0,zeroVector);

    // Make the poly-liner
    G4Box* boxMagnetPLine = new G4Box("boxMagnetPLine", stmMagnetHoleHalfLengths[0]+pSTMShieldPipeParams.linerWidth(), stmMagnetHoleHalfLengths[1]+pSTMShieldPipeParams.linerWidth(), pSTMMagnetParams.zHalfLength());//wyq change the length of liner,also next line need to make longer
    G4Box* boxPolyHole = new G4Box("boxPolyHole", stmMagnetHoleHalfLengths[0], stmMagnetHoleHalfLengths[1], stmMagnetHoleHalfLengths[2]+1.0);

    VolumeInfo stmMagnetPLine;
    stmMagnetPLine.name = "stmMagnetPLine";
    stmMagnetPLine.solid = new G4SubtractionSolid(stmMagnetPLine.name,boxMagnetPLine,boxPolyHole,0,zeroVector);



    if (pSTMMagnetParams.build()){
      finishNesting(stmMagnet,
                    findMaterialOrThrow(pSTMMagnetParams.materialName()),
                    0,
                    stmMagnetPositionInParent, //mstmMagnetPositionInMother,
                    parentInfo.logical, //parentInfo.logical, //mstmMotherInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Gray(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);

      if(pSTMMagnetParams.materialName()!="G4_AIR")//by Yaqian for removing the magnet and poly
			finishNesting(stmMagnetPLine,
                    findMaterialOrThrow(pSTMShieldPipeParams.materialLiner()),
                    0,
                    stmMagnetPositionInParent, //mstmMagnetPositionInMother,
                    parentInfo.logical, //parentInfo.logical, //mstmMotherInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Gray(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);


    }

    if ( verbosityLevel > 0) {
       cout << __func__ << " Sweeper magnet extent in z   : "
            << pSTMMagnetParams.zBegin() <<","<< pSTMMagnetParams.zEnd() << endl;
       cout << __func__ << " Sweeper magnet hole opening xHalfLength : "<< stmMagnetHoleHalfLengths[0] << endl;
       cout << __func__ << " Sweeper magnet hole opening yHalfLength : "<< stmMagnetHoleHalfLengths[1] << endl;
    }


    //===================== Transport Pipe ==========================

    //create a disk so we can subtract a space for the existing VD to fit inside (avoid overlaps)
    GeomHandle<VirtualDetector> vdg;
    const double vdHL = CLHEP::mm * vdg->getHalfLength();
    const double vdR  = _config.getDouble("vd.DSNeutronShieldExit.r");
    G4Tubs* aDiskVDDSNeutronShieldExitTub = new G4Tubs( "subtSpaceForVDDSNeutronShieldExit",
                                       0.0,
                                       vdR+0.01,
                                       vdHL+0.01,// a bit larger to avoid overlap with VD
                                       0.0,
                                       CLHEP::twopi);
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionInMu2e = vdg->getGlobal(VirtualDetectorId::DSNeutronShieldExit);
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionInParent = vdDSNeutronShieldExitPositionInMu2e - parentCenterInMu2e;
    if (verbosityLevel>0){
      std::cout << __func__ << " vdDSNeutronShieldExitPositionInMu2e = "<<vdDSNeutronShieldExitPositionInMu2e<<std::endl;
    }


    //create a disk so we can subtract a space for the existing VD to fit inside (avoid overlaps)
    G4Tubs* aDiskVDSTM_UpStrTub = new G4Tubs( "subtSpaceForVDSTM_UpStr",
                                       0.0,
                                       1000.0, //make is very large to be sure we subtract enough
                                       vdHL+0.01,// a bit larger to avoid overlap with VD
                                       0.0,
                                       CLHEP::twopi);
    CLHEP::Hep3Vector vdSTM_UpStrPositionInMu2e = vdg->getGlobal(VirtualDetectorId::STM_UpStr);
    CLHEP::Hep3Vector vdSTM_UpStrPositionInParent = vdSTM_UpStrPositionInMu2e - parentCenterInMu2e;
    if (verbosityLevel>0){
      std::cout << __func__ << " vdSTM_UpStrPositionInMu2e = "<<vdSTM_UpStrPositionInMu2e<<std::endl;
    }

    //create a disk so we can subtract a space for the existing VD to fit inside (avoid overlaps)
//     G4Tubs* aDiskVDSTM_CRVShieldDnStrTub = new G4Tubs( "subtSpaceForVDSTM_CRVShieldDnStr",
//                                        0.0,
//                                        1000.0, //make is very large to be sure we subtract enough
//                                        vdHL+0.01,// a bit larger to avoid overlap with VD
//                                        0.0,
//                                        CLHEP::twopi);
//     CLHEP::Hep3Vector vdSTM_CRVShieldDnStrPositionInMu2e   = vdg->getGlobal(VirtualDetectorId::STM_CRVShieldDnStr);
//     CLHEP::Hep3Vector vdSTM_CRVShieldDnStrPositionInParent = vdSTM_CRVShieldDnStrPositionInMu2e - parentCenterInMu2e;
//     if (verbosityLevel>0){
//       std::cout<<"vdSTM_CRVShieldDnStrPositionInMu2e = "<<vdSTM_CRVShieldDnStrPositionInMu2e<<std::endl;
//     }


    //---


    const double flangeHalfLength      =     pSTMTransportPipeParams.flangeHalfLength();
    const double flangeFullLength      = 2.0*pSTMTransportPipeParams.flangeHalfLength();

    //make the region of pipe that will contain the magnetic field
    G4Tubs* aPipeCenterTub = new G4Tubs("aPipeCenterTub",
					pSTMTransportPipeParams.radiusIn(),  //inner radius
					pSTMTransportPipeParams.radiusOut(), //outer radius
					stmMagnetHoleHalfLengths[2],
					0.0,
					CLHEP::twopi);

    VolumeInfo pipeCenterTubInfo;
    pipeCenterTubInfo.name = "pipeCenterTub";
    pipeCenterTubInfo.solid = aPipeCenterTub;

    //make the gas that goes inside the region of pipe that contains the magnetic field
    G4Tubs* aPipeCenterGasTub = new G4Tubs("aPipeCenterGasTub",
                                           0.0, //inner radius of gas
                                           pSTMTransportPipeParams.radiusIn(), //outer radius of gas
                                           stmMagnetHoleHalfLengths[2],
                                           0.0,
                                           CLHEP::twopi);
    VolumeInfo pipeCenterGasTubInfo;
    pipeCenterGasTubInfo.name = "pipeGasTub";
    pipeCenterGasTubInfo.solid = aPipeCenterGasTub;

    //make the region of pipe that goes downstream of the magnetic field, with a flange
    G4Tubs* aPipeDnStrTub = new G4Tubs("aPipeDnStrTub",
                                       pSTMTransportPipeParams.radiusIn(),  //inner radius
                                       pSTMTransportPipeParams.radiusOut()+pSTMTransportPipeParams.flangeOverhangR(), //outer radius
                                       pSTMTransportPipeParams.dnStrHalflength(),
                                       0.0,
                                       CLHEP::twopi);
    G4Tubs* aPipeDnStrSubtTub = new G4Tubs("aPipeDnStrSubtTub",
                                           pSTMTransportPipeParams.radiusOut(), //inner radius to subtract
                                           pSTMTransportPipeParams.radiusOut()+2.0*pSTMTransportPipeParams.flangeOverhangR(), //outer radius, make sure to subtract enough
                                           pSTMTransportPipeParams.dnStrHalflength(),
                                           0.0,
                                           CLHEP::twopi);
    CLHEP::Hep3Vector flangeOffset(0.0, 0.0, -flangeFullLength);
    VolumeInfo pipeDnStrTubInfo;
    pipeDnStrTubInfo.name = "pipeDnStrTub";
    pipeDnStrTubInfo.solid = new G4SubtractionSolid(pipeDnStrTubInfo.name,aPipeDnStrTub, aPipeDnStrSubtTub, 0, flangeOffset);
    CLHEP::Hep3Vector pipeDnStrOffset(0.0, 0.0, stmMagnetHoleHalfLengths[2]+pSTMTransportPipeParams.dnStrHalflength());

    //make a downstream window to hold the helium or vacuum inside the pipe
    G4Tubs* aPipeDnStrWindowTub = new G4Tubs( "aPipeDnStrWindowTub",
                                              0.0, // inner radius of window
                                              pSTMTransportPipeParams.radiusIn(), //outer radius of window
                                              pSTMTransportPipeParams.dnStrWindowHalflength(),
                                              0.0,
                                              CLHEP::twopi);
    VolumeInfo pipeDnStrWindowTubInfo;
    pipeDnStrWindowTubInfo.name = "pipeDnStrWindowTub";
    pipeDnStrWindowTubInfo.solid = aPipeDnStrWindowTub;
    CLHEP::Hep3Vector pipeDnStrWindowOffset(0.0, 0.0, stmMagnetHoleHalfLengths[2]+2.0*pSTMTransportPipeParams.dnStrHalflength()-flangeHalfLength);

    //put gas/vacuum inside the downstream portion of the pipe
    const double gasDnStrHalfLength = 0.5*(2.0*pSTMTransportPipeParams.dnStrHalflength() - flangeHalfLength - pSTMTransportPipeParams.dnStrWindowHalflength() );
    G4Tubs* aPipeGasDnStrTub = new G4Tubs( "aPipeGasDnStrTub",
                                            0.0, //inner radius
                                            pSTMTransportPipeParams.radiusIn(), //outer radius of gas
                                            gasDnStrHalfLength,
                                            0.0,
                                            CLHEP::twopi);
    VolumeInfo pipeGasDnStrTubInfo;
    pipeGasDnStrTubInfo.name = "pipeGasDnStrTub";
    pipeGasDnStrTubInfo.solid = aPipeGasDnStrTub;
    CLHEP::Hep3Vector pipeGasDnStrOffset(0.0, 0.0, stmMagnetHoleHalfLengths[2]+gasDnStrHalfLength);

    //make the region of pipe that goes upstream of the magnetic field, with a flange
    const double IFB_endplug_z_center     = ds->cryoZMax() + _config.getDouble("ifb.endplug.z");
    const double IFB_endplug_z_halflength = _config.getDouble("ifb.endplug.halfLength");
    const double pipeUpStrHalfLength = 0.5*((stmMagnetPositionInMu2e.z()-stmMagnetHoleHalfLengths[2]) - (IFB_endplug_z_center+IFB_endplug_z_halflength)) - pSTMTransportPipeParams.upStrSpace();//leave a space between pipe and IFB

    G4Tubs* aPipeUpStrTub     = new G4Tubs("aPipeUpStrTub",
                                           pSTMTransportPipeParams.radiusIn(), //inner radius
                                           pSTMTransportPipeParams.radiusOut()+pSTMTransportPipeParams.flangeOverhangR(), //outer radius
                                           pipeUpStrHalfLength,//
                                           0.0,//
                                           CLHEP::twopi);
    G4Tubs* aPipeUpStrSubtTub = new G4Tubs("aPipeUpStrSubtTub",
                                           pSTMTransportPipeParams.radiusOut(), //inner radius
                                           pSTMTransportPipeParams.radiusOut()+2.0*pSTMTransportPipeParams.flangeOverhangR(), //outer radius, make sure to subtract enough
                                           pipeUpStrHalfLength,//
                                           0.0,//
                                           CLHEP::twopi);
    CLHEP::Hep3Vector pipeUpStrOffset(0.0, 0.0, -stmMagnetHoleHalfLengths[2]-pipeUpStrHalfLength);
    //first make the subtraction so it has a flange
    G4SubtractionSolid *pipeUpStrTubTemp1 = new G4SubtractionSolid("pipeUpStrTubTemp1",aPipeUpStrTub, aPipeUpStrSubtTub, 0, -flangeOffset);

    //subtract a slice so VDDSNeutronShieldExit can fit through the pipe without overlap
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionWRTpipeUpStr = vdDSNeutronShieldExitPositionInMu2e - (stmMagnetPositionInMu2e+pipeUpStrOffset);
    CLHEP::Hep3Vector vdSTM_UpStrPositionWRTpipeUpStr = vdSTM_UpStrPositionInMu2e - (stmMagnetPositionInMu2e+pipeUpStrOffset);
    //std::cout<<"vdDSNeutronShieldExitPositionWRTpipe = "<<vdDSNeutronShieldExitPositionWRTpipeUpStr<<std::endl;
    //std::cout<<"vdSTM_UpStrPositionWRTpipe = "<<vdDSNeutronShieldExitPositionWRTpipeUpStr<<std::endl;
    G4SubtractionSolid *pipeUpStrTubTemp2 = new G4SubtractionSolid("pipeUpStrTubTemp2",pipeUpStrTubTemp1, aDiskVDDSNeutronShieldExitTub,      0, vdDSNeutronShieldExitPositionWRTpipeUpStr);
    G4SubtractionSolid *pipeUpStrTubTemp3 = new G4SubtractionSolid("pipeUpStrTubTemp3",pipeUpStrTubTemp2, aDiskVDSTM_UpStrTub,      0, vdSTM_UpStrPositionWRTpipeUpStr);
		//>Yaqian Wang insert Poly Up and Dnstream of IFB Window
		double rSTMIFBLiner=_config.getDouble("ifb.endwindow.rOut");
    // G4Tubs *tubIFBLinerUp = new G4Tubs("tubIFBLinerUp", 0.0, rSTMIFBLiner, _config.getDouble("stm.FOVcollimator.absorber.halfLength"), 0.0, CLHEP::twopi );
    // G4Tubs *tubIFBLinerDn = new G4Tubs("tubIFBLinerDn", 0.0, rSTMIFBLiner, _config.getDouble("stm.FOVcollimator.absorber.halfLength"), 0.0, CLHEP::twopi );
    G4Tubs *tubIFBLinerUp = new G4Tubs("tubIFBLinerUp", 0.0, rSTMIFBLiner, _config.getDouble("stm.IFBWindowLinerUp.halfLength"), 0.0, CLHEP::twopi );
    G4Tubs *tubIFBLinerDn = new G4Tubs("tubIFBLinerDn", 0.0, rSTMIFBLiner, _config.getDouble("stm.IFBWindowLinerDn.halfLength"), 0.0, CLHEP::twopi );
    VolumeInfo STMIFBLinerUp;
    STMIFBLinerUp.name="STMIFBLinerUp";
    STMIFBLinerUp.solid=tubIFBLinerUp;
    VolumeInfo STMIFBLinerDn;
    STMIFBLinerDn.name="STMIFBLinerDn";
    STMIFBLinerDn.solid=tubIFBLinerDn;
	VolumeInfo const & ifbParentInfo    = _helper->locateVolInfo("DS3Vacuum");
    // VolumeInfo const & ifbEndWindowInfo = _helper->locateVolInfo("IFB_EndWindow");
    VolumeInfo const & ifbEndPlug       = _helper->locateVolInfo("IFB_EndPlug");
    G4ThreeVector STMIFBLinerUpPositionInParent = G4ThreeVector(0,0,17375-26);
    G4ThreeVector STMIFBLinerDnPositionInParent = G4ThreeVector(0,0,17415-26);
//    G4ThreeVector STMIFBLinerDnPositionInParent = pSTMMagnetParams.originInMu2e() -G4ThreeVector(0,0,pSTMMagnetParams.zHalfLength()+2*pSTMShieldPipeParams.pipeHalfLength()+pSTMShieldPipeParams.dnStrWallHalflength()*2+1150);
//		const VolumeInfo& mbsParent = _helper->locateVolInfo("MBSMother");

    // Downstream liner needs to be in other volume...


if(_config.getBool("stm.IFBWindowLinerUp.build",false)){
    // G4ThreeVector posInParent =
    //   ifbEndWindowInfo.centerInParent -
    //   G4ThreeVector(0.,
    //                 0.,
    //                 _config.getDouble("stm.IFBWindowLinerUp.halfLength") +
    //                 _config.getDouble("ifb.endwindow.halfLength"));
    finishNesting(STMIFBLinerUp,
       findMaterialOrThrow(_config.getString("stm.IFBWindowLinerUp.material")),
       0x0,
       STMIFBLinerUpPositionInParent,
       // posInParent,
       ifbParentInfo.logical,
       0,
       STMisVisible,
       G4Color::Gray(),
       STMisSolid,
       forceAuxEdgeVisible,
       placePV,
       doSurfaceCheck
       );
}
if(_config.getBool("stm.IFBWindowLinerDn.build",false)){
    // G4ThreeVector posInParent =
    //   ifbEndWindowInfo.centerInParent +
    //   G4ThreeVector(0.,
    //                 0.,
    //                 _config.getDouble("stm.IFBWindowLinerDn.halfLength") +
    //                 _config.getDouble("ifb.endwindow.halfLength"));
    G4ThreeVector EndPlugPosInParent = ifbEndPlug.centerInMu2e()-ifbParentInfo.centerInMu2e();
    G4ThreeVector posInParent = EndPlugPosInParent +
      G4ThreeVector(0.,
                    0.,
                    _config.getDouble("stm.IFBWindowLinerDn.halfLength") +
                    _config.getDouble("ifb.endplug.halfLength")) +
      G4ThreeVector(0., 0., _config.getDouble("stm.IFBWindowLinerDn.zshift"));
    // G4ThreeVector posInParent =
    //   ifbParentInfo.centerInParent +
    //   G4ThreeVector(0.,
    //                 0.,
    //                 _config.getDouble("stm.IFBWindowLinerDn.halfLength") +
    //                 _config.getDouble("ifb.endplug.halfLength"));
    STMIFBLinerDnPositionInParent +=
      G4ThreeVector(0., 0., _config.getDouble("stm.IFBWindowLinerDn.zshift"));
    // G4ThreeVector posInParent =
    //   ifbEndWindowInfo.centerInParent +
    //   G4ThreeVector(0.,
    //                 0.,
    //                 500.);
    finishNesting(STMIFBLinerDn,
       findMaterialOrThrow(_config.getString("stm.IFBWindowLinerDn.material")),
       0x0,
       STMIFBLinerDnPositionInParent,
       // posInParent,
       // parentInfo.logical,
       ifbParentInfo.logical,
       0,
       STMisVisible,
       G4Color::Red(),
       STMisSolid,
       forceAuxEdgeVisible,
       placePV,
       doSurfaceCheck
       );
}



	//<
	 //>Yaqian Wang insert poly in the concrete ext hole
		double ExtShieldlinerHalfLength=_config.getDouble("stm.ExtShieldLiner.HalfLength");
		double ExtShieldlinerRIn=_config.getDouble("stm.ExtShieldLiner.RIn");
//w		double ExtShieldlinerROut=_config.getDouble("stm.ExtShieldLiner.ROut");
		double ExtShieldlinerROut=_config.getDouble("ExtShieldDownstream.holeRadiusType11Box4Hole1");

    G4Tubs* ExtShieldLiner=new G4Tubs("ExtShieldLiner",
			 ExtShieldlinerRIn, //inner radius
			 ExtShieldlinerROut,//outer radius
			 ExtShieldlinerHalfLength,//half length
			 0.0,//
			 CLHEP::twopi);
    const double pipeCrack = _config.getDouble("stm.shield.pipe.crack",0.0);
		VolumeInfo ExtShieldLinerInfo;
		ExtShieldLinerInfo.name="ExtShieldLiner";
		ExtShieldLinerInfo.solid=ExtShieldLiner;
    G4ThreeVector stmExtShieldLinerPositionInParent = pSTMMagnetParams.originInMu2e() - parentCenterInMu2e-G4ThreeVector(0,0,pSTMMagnetParams.zHalfLength()+2*pSTMShieldPipeParams.pipeHalfLength()+pSTMShieldPipeParams.dnStrWallHalflength()*2+ExtShieldlinerHalfLength);//5 is crack dn concrete
		vector<double> tolsType11;//tolerance, 2nd is for y,rotate later
		_config.getVectorDouble("ExtShieldDownstream.tolsType11",tolsType11,3);
		vector<double> outlineType11VVerts;
		_config.getVectorDouble("ExtShieldDownstream.outlineType11VVerts",outlineType11VVerts,8);//1st member is -454, half of the length
		double concreteHalfLength=-1.0*outlineType11VVerts[0]+0.5*tolsType11[1];
    G4Tubs* ExtShieldLinerC=new G4Tubs("ExtShieldLinerC",
			 ExtShieldlinerRIn, //inner radius
			 ExtShieldlinerROut,//outer radius
			 concreteHalfLength-ExtShieldlinerHalfLength,//half length
			 0.0,//
			 CLHEP::twopi);
		VolumeInfo ExtShieldLinerCInfo;
		ExtShieldLinerCInfo.name="ExtShieldLinerC";
		ExtShieldLinerCInfo.solid=ExtShieldLinerC;
    G4ThreeVector stmExtShieldLinerCPositionInParent =stmExtShieldLinerPositionInParent-G4ThreeVector(0,0,concreteHalfLength);

if(_config.getBool("stm.ExtShieldLiner.build")){
    finishNesting(ExtShieldLinerInfo,
       findMaterialOrThrow(pSTMShieldPipeParams.materialLiner()),
       0x0,
       stmExtShieldLinerPositionInParent,
       parentInfo.logical,
       0,
       STMisVisible,
       G4Color::Gray(),
       STMisSolid,
       forceAuxEdgeVisible,
       placePV,
       doSurfaceCheck
       );
     finishNesting(ExtShieldLinerCInfo,//fill concrete back
        findMaterialOrThrow("CONCRETE_MARS"),
        0x0,
        stmExtShieldLinerCPositionInParent,
        parentInfo.logical,
        0,
        STMisVisible,
        G4Color::Gray(),
        STMisSolid,
        forceAuxEdgeVisible,
        placePV,
        doSurfaceCheck
        );
}
	//<
    VolumeInfo pipeUpStrTubInfo;
    pipeUpStrTubInfo.name = "pipeUpStrTub";
    pipeUpStrTubInfo.solid = pipeUpStrTubTemp3;


    //put gas inside the upstream portion of the pipe
    const double pipeGasUpStrHalfLength = 0.5*(2.0*pipeUpStrHalfLength - flangeHalfLength - pSTMTransportPipeParams.upStrWindowHalflength());
    G4Tubs* aPipeGasUpStrTub = new G4Tubs( "aPipeGasUpStrTub",
                                            0.0, //inner radius
                                            pSTMTransportPipeParams.radiusIn(), //outer radius
                                            pipeGasUpStrHalfLength,//
                                            0.0,//
                                            CLHEP::twopi);
    CLHEP::Hep3Vector pipeGasUpStrOffset(0.0, 0.0, -stmMagnetHoleHalfLengths[2]-pipeGasUpStrHalfLength);
    //subtract a slice so VDDSNeutronShieldExit can fit through the pipe without overlap
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionWRTpipeGasUpStr = vdDSNeutronShieldExitPositionInMu2e - (stmMagnetPositionInMu2e+pipeGasUpStrOffset);
    CLHEP::Hep3Vector vdSTM_UpStrPositionWRTpipeGasUpStr = vdSTM_UpStrPositionInMu2e - (stmMagnetPositionInMu2e+pipeGasUpStrOffset);
    //std::cout<<"vdDSNeutronShieldExitPositionWRTtube = "<<vdDSNeutronShieldExitPositionWRTpipeGasUpStr<<std::endl;
    //std::cout<<"vdSTM_UpStrPositionWRTtube = "<<vdSTM_UpStrPositionWRTpipeGasUpStr<<std::endl;
    G4SubtractionSolid *pipeGasUpStrTubTemp1 = new G4SubtractionSolid("pipeGasUpStrTubTemp1",aPipeGasUpStrTub,     aDiskVDDSNeutronShieldExitTub, 0, vdDSNeutronShieldExitPositionWRTpipeGasUpStr);
    G4SubtractionSolid *pipeGasUpStrTubTemp2 = new G4SubtractionSolid("pipeGasUpStrTubTemp2",pipeGasUpStrTubTemp1, aDiskVDSTM_UpStrTub, 0, vdSTM_UpStrPositionWRTpipeGasUpStr);
    VolumeInfo pipeGasUpStrTubInfo;
    pipeGasUpStrTubInfo.name = "pipeGasUpStrTub";
    pipeGasUpStrTubInfo.solid = pipeGasUpStrTubTemp2;

    //make an upstream window to hold the gas/vacuum in the pipe
    G4Tubs* aPipeUpStrWindowTub = new G4Tubs( "aPipeUpStrWindowTub",
                                            0.0, // inner radius
                                            pSTMTransportPipeParams.radiusIn(), //outer radius of window
                                            pSTMTransportPipeParams.upStrWindowHalflength(),
                                            0.0,
                                            CLHEP::twopi);
    VolumeInfo pipeUpStrWindowTubInfo;
    pipeUpStrWindowTubInfo.name = "pipeUpStrWindowTub";
    pipeUpStrWindowTubInfo.solid = aPipeUpStrWindowTub;
    CLHEP::Hep3Vector pipeUpStrWindowOffset(0.0, 0.0, -stmMagnetHoleHalfLengths[2]-2.0*pipeUpStrHalfLength+flangeHalfLength);


    if (pSTMTransportPipeParams.build()){
      if (verbosityLevel>0){
        const double pipeTotalHalfLength = pipeUpStrHalfLength+stmMagnetHoleHalfLengths[2]+pSTMTransportPipeParams.dnStrHalflength();
        std::cout<<__func__<<" STM Transport Pipe z_halflength = "<< pipeTotalHalfLength <<std::endl;
        std::cout<<__func__<<" STM Transport Pipe z_min        = "<< (stmMagnetPositionInMu2e.z()+pipeUpStrOffset.z())-pipeUpStrHalfLength <<std::endl;
        std::cout<<__func__<<" STM Transport Pipe z_max        = "<< (stmMagnetPositionInMu2e.z()+pipeDnStrOffset.z())+pSTMTransportPipeParams.dnStrHalflength() <<std::endl;
        std::cout<<__func__<<" STM Transport Pipe rIn          = "<< pSTMTransportPipeParams.radiusIn() <<std::endl;
      }
      finishNesting(pipeCenterTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.material()),
                0x0,
                stmMagnetPositionInParent,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeCenterGasTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.gasMaterial()),
                0x0,
                stmMagnetPositionInParent,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeDnStrTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.material()),
                0x0,
                stmMagnetPositionInParent+pipeDnStrOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeDnStrWindowTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.dnStrWindowMaterial()),
                0x0,
                stmMagnetPositionInParent+pipeDnStrWindowOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeGasDnStrTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.gasMaterial()),
                0x0,
                stmMagnetPositionInParent+pipeGasDnStrOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeUpStrTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.material()),
                0x0,
                stmMagnetPositionInParent+pipeUpStrOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeGasUpStrTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.gasMaterial()),
                0x0,
                stmMagnetPositionInParent+pipeGasUpStrOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
      finishNesting(pipeUpStrWindowTubInfo,
                findMaterialOrThrow(pSTMTransportPipeParams.upStrWindowMaterial()),
                0x0,
                stmMagnetPositionInParent+pipeUpStrWindowOffset,
                parentInfo.logical,
                0,
                STMisVisible,
                G4Color::Blue(),
                STMisSolid,
                forceAuxEdgeVisible,
                placePV,
                doSurfaceCheck
                );
    }



    // ----- Magnetic Field -----------------------------------------------------
    //Create a magnetic field inside the window (hole) of the magnet box
    //and in the pipe, and pipe gas, that goes through the magnet
    //Note the local values for the stepper etc...
    //Geant4 should take ownership of the objects created here

    const double stmMagnetFieldHalfLengths[3] = {pSTMMagnetParams.xHoleHalfLength(),
                                                pSTMMagnetParams.yHoleHalfLength(),
                                                pSTMMagnetParams.zHalfLength()};
//wyq remove crack                               pSTMMagnetParams.zHalfLength()-pSTMShieldPipeParams.linerWidth()};

    VolumeInfo stmMagneticFieldBoxInfo;
    stmMagneticFieldBoxInfo.name = "stmMagneticField";

    // Make another rectangular volume for the magnetic field
    G4Box* boxField  = new G4Box("boxField",stmMagnetFieldHalfLengths[0],stmMagnetFieldHalfLengths[1],stmMagnetFieldHalfLengths[2]);

    G4Tubs* aPipeCenterTubSubt = new G4Tubs( "aPipeCenterTubSubt",
                                 0.0, //inner radius 0.0 to subtract also the gas region
                                 pSTMTransportPipeParams.radiusOut()+0.01, //outer radius
                                 stmMagnetFieldHalfLengths[2]+1.0,// make the subtraction slightly larger to avoid edge effects
                                 0.0,
                                 CLHEP::twopi);

    if (pSTMMagnetParams.build()){
      if  (pSTMTransportPipeParams.build()){

          stmMagneticFieldBoxInfo.solid = new G4SubtractionSolid(stmMagneticFieldBoxInfo.name,boxField,aPipeCenterTubSubt,0,zeroVector);
          finishNesting(stmMagneticFieldBoxInfo,
                        findMaterialOrThrow(_config.getString("hall.insideMaterialName")),
                        0x0,
                        stmMagnetPositionInParent,
                        parentInfo.logical,
                        0,
                        pSTMMagnetParams.fieldVisible(),   //magnetic field itself visible
                        G4Color::Blue(),
                        false,
                        forceAuxEdgeVisible,
                        placePV,
                        doSurfaceCheck
                       );

      } else {
          stmMagneticFieldBoxInfo = nestBox("stmMagneticField",
                                            stmMagnetFieldHalfLengths,
                                            findMaterialOrThrow(_config.getString("hall.insideMaterialName")),
                                            0x0,
                                            stmMagnetPositionInParent,
                                            parentInfo,
                                            0,
                                            pSTMMagnetParams.fieldVisible(),
                                            G4Color::Blue(),
                                            false,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                           );
      }

      //G4MagneticField        *localMagField        = new G4UniformMagField(G4ThreeVector(pSTMMagnetParams.field()*CLHEP::tesla,0.0,0.0));//This makes negatively charged particles go towards the floor
      G4MagneticField        *localMagField        = new G4UniformMagField(G4ThreeVector(pSTMMagnetParams.fieldx()*CLHEP::tesla,pSTMMagnetParams.fieldy()*CLHEP::tesla,pSTMMagnetParams.fieldz()*CLHEP::tesla));
      G4Mag_EqRhs            *MagRHS               = new G4Mag_UsualEqRhs(localMagField);
      G4MagIntegratorStepper *localMagStepper      = new G4ExactHelixStepper(MagRHS); // we use a specialized stepper
      G4ChordFinder          *localMagChordFinder  = new G4ChordFinder(localMagField,1.0e-2*CLHEP::mm,localMagStepper);
      G4FieldManager         *localMagFieldManager = new G4FieldManager(localMagField,localMagChordFinder,false);// pure magnetic filed does not change energy

      stmMagneticFieldBoxInfo.logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
      if (pSTMTransportPipeParams.build()){
        pipeCenterTubInfo.logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
        pipeCenterGasTubInfo.logical->SetFieldManager(localMagFieldManager, true); // last "true" arg propagates field to all volumes it contains
      }
      G4UserLimits* mstmMagStepLimit = new G4UserLimits(5.*CLHEP::mm);
      stmMagneticFieldBoxInfo.logical->SetUserLimits(mstmMagStepLimit);
      if (pSTMTransportPipeParams.build()){
        pipeCenterTubInfo.logical->SetUserLimits(mstmMagStepLimit);
        pipeCenterGasTubInfo.logical->SetUserLimits(mstmMagStepLimit);
      }

    }


    //===================== Field-of-View (FOV) Collimator ==========================

    const double stmFOVCollHalfLength1 = pSTMFOVCollimatorParams.halfLength();
    const double stmFOVCollHalfWidth1  = pSTMFOVCollimatorParams.halfWidth();
    const double stmFOVCollHalfHeight1 = pSTMFOVCollimatorParams.halfHeight();
    const double stmFOVCollHalfLength2 = pSTMFOVCollimatorParams.linerHalfLength();
    const double stmFOVCollHalfWidth2  = pSTMFOVCollimatorParams.linerHalfWidth();
    const double stmFOVCollHalfHeight2 = pSTMFOVCollimatorParams.linerHalfHeight();

    // position of FOV collimator
    G4ThreeVector stmFOVCollPositionInMu2e1   = pSTMFOVCollimatorParams.originInMu2e();
    G4ThreeVector stmFOVCollPositionInParent1 = pSTMFOVCollimatorParams.originInMu2e() - parentCenterInMu2e;
    // Make the box for the collimator
    G4Box* boxFOVColl = new G4Box("boxFOVColl",stmFOVCollHalfWidth1,stmFOVCollHalfHeight1,stmFOVCollHalfLength1);
    //Make the tube for the hole
    G4Tubs *tubFOVColl1 = new G4Tubs("tubFOVColl1", 0.0, pSTMFOVCollimatorParams.hole1RadiusUpStr(), stmFOVCollHalfLength1+1.0, 0.0, CLHEP::twopi );
    //Make a box to subtract so liner can fit inside
    G4Box* boxFOVCollLinerToSubt = new G4Box("boxFOVCollLinerToSubt",stmFOVCollHalfWidth2+1.0,stmFOVCollHalfHeight2+1.0,pSTMFOVCollimatorParams.linerCutOutHalfLength()+1.0);
    // Combine into the collimator with the liner cutout and collimation hole
    VolumeInfo collimatorFOV;
    collimatorFOV.name = "collimatorFOV";
    collimatorFOV.solid = new G4SubtractionSolid(collimatorFOV.name,boxFOVColl,         tubFOVColl1,0,G4ThreeVector(0.0,0.0,0.0));
    collimatorFOV.solid = new G4SubtractionSolid(collimatorFOV.name,collimatorFOV.solid,boxFOVCollLinerToSubt,0,G4ThreeVector(0.0,0.0,-1.0*(stmFOVCollHalfLength1-pSTMFOVCollimatorParams.linerCutOutHalfLength() )));

    //position of liner
    G4ThreeVector stmFOVCollPositionInMu2e2   = stmFOVCollPositionInMu2e1   + G4ThreeVector(0.0,0.0, -stmFOVCollHalfLength1+2.0*pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);
    G4ThreeVector stmFOVCollPositionInParent2 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,0.0, -stmFOVCollHalfLength1+2.0*pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);
    // make the box for the liner
    G4Box* boxFOVCollLiner = new G4Box("boxFOVCollLiner",stmFOVCollHalfWidth2,stmFOVCollHalfHeight2,stmFOVCollHalfLength2);
    //Make the tube for the hole
    G4Tubs *tubFOVCollLiner1 = new G4Tubs("tubFOVCollLiner1", 0.0, pSTMFOVCollimatorParams.hole1RadiusUpStr(), stmFOVCollHalfLength2+1.0, 0.0, CLHEP::twopi );
    // Combine into the liner with hole
    VolumeInfo collimatorFOVliner;
    collimatorFOVliner.name = "collimatorFOVliner";
    collimatorFOVliner.solid = new G4SubtractionSolid(collimatorFOVliner.name,boxFOVCollLiner,tubFOVCollLiner1,0,G4ThreeVector(0.0,0.0,0.0));

    // Liner sheets to cover the upstream corner of FOV collimator
    VolumeInfo FOVlinerH;
    FOVlinerH.name = "FOVlinerH";
    FOVlinerH.solid = new G4Box(FOVlinerH.name,stmFOVCollHalfWidth2, pSTMShieldPipeParams.linerWidth()/2.0, pSTMFOVCollimatorParams.linerCutOutHalfLength()-stmFOVCollHalfLength2);
    //G4ThreeVector stmFOVCollPositionInParent3 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,-stmFOVCollHalfHeight2+pSTMShieldPipeParams.linerWidth()/2.0, -2*stmFOVCollHalfLength2);
    G4ThreeVector stmFOVCollPositionInParent3 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,-stmFOVCollHalfHeight2+pSTMShieldPipeParams.linerWidth()/2.0, -2*stmFOVCollHalfLength2);
/*remove by Yaqian Since no crack between magnet and FOV, thus no space for this
    VolumeInfo FOVlinerW;
    FOVlinerW.name = "FOVlinerW";
    double cornerHeight = (stmFOVCollHalfLength1-stmFOVCollHalfLength2)/2.0;
		double FOVlinerWHalfLength=pSTMShieldPipeParams.linerWidth()/2.0;
    FOVlinerW.solid = new G4Box(FOVlinerW.name,stmFOVCollHalfWidth2, cornerHeight/2.0, FOVlinerWHalfLength);
		G4ThreeVector stmFOVCollPositionInParent4 = stmFOVCollPositionInParent1 + G4ThreeVector(0.0,-stmFOVCollHalfHeight2-cornerHeight/2.0, -stmFOVCollHalfLength1-FOVlinerWHalfLength);
		*/
		double FOVlinerW2HalfLength=_config.getDouble("stm.FOVcollimator.linerWHalfLength");//in Magnetic field
    VolumeInfo const & STMagInfo = _helper->locateVolInfo("stmMagneticField");
		if(FOVlinerW2HalfLength!=0){
    VolumeInfo FOVlinerW2;
    FOVlinerW2.name = "FOVlinerW2";
    //double cornerHeight = (stmFOVCollHalfLength1-stmFOVCollHalfLength2)/2.0;
    double cornerHeight = stmFOVCollHalfHeight1-stmFOVCollHalfHeight2-(pSTMMagnetParams.yHalfLength()-stmMagnetHoleHalfLengths[1])+2*pSTMShieldPipeParams.linerWidth();//the holelength include 1cm liner
    FOVlinerW2.solid=new G4Box(FOVlinerW2.name,pSTMMagnetParams.xHoleHalfLength()-0.01,cornerHeight/2.0,FOVlinerW2HalfLength);
		G4ThreeVector stmFOVCollPositionInParent5=G4ThreeVector(0.0,-pSTMMagnetParams.yHoleHalfLength()+(cornerHeight)/2.0,pSTMMagnetParams.zHalfLength()-FOVlinerW2HalfLength);
      finishNesting(FOVlinerW2,
                    findMaterialOrThrow(pSTMFOVCollimatorParams.linerMaterial()),
                    0,
                    stmFOVCollPositionInParent5,
                    STMagInfo.logical,
                    0,
                    true,
                    G4Colour::Magenta(),
                    false,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
}
    if (pSTMFOVCollimatorParams.build()){
      finishNesting(collimatorFOV,
                    findMaterialOrThrow(pSTMFOVCollimatorParams.material()),
                    0,
                    stmFOVCollPositionInParent1,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
      finishNesting(collimatorFOVliner,
                    findMaterialOrThrow(pSTMFOVCollimatorParams.linerMaterial()),
                    0,
                    stmFOVCollPositionInParent2,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
      finishNesting(FOVlinerH,
                    findMaterialOrThrow(pSTMFOVCollimatorParams.linerMaterial()),
                    0,
                    stmFOVCollPositionInParent3,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
/*remove by Yaqian see the above lines where the solid is made
      finishNesting(FOVlinerW,
                    findMaterialOrThrow(pSTMFOVCollimatorParams.linerMaterial()),
                    0,
                    stmFOVCollPositionInParent4,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
										*/

    }

    G4Tubs *tubFOVCollAbsorber = new G4Tubs("tubFOVCollAbsorber", 0.0, pSTMFOVCollimatorParams.hole1RadiusUpStr()-0.01, _config.getDouble("stm.FOVcollimator.absorber.halfLength"), 0.0, CLHEP::twopi );
    VolumeInfo collimatorFOVAbsorber;
    collimatorFOVAbsorber.name = "collimatorFOVAbsorber";
    collimatorFOVAbsorber.solid = tubFOVCollAbsorber;
    G4ThreeVector stmFOVCollAbsorberPositionInParent = stmFOVCollPositionInParent2 + G4ThreeVector(0.0,0.0, stmFOVCollHalfLength2-_config.getDouble("stm.FOVcollimator.absorber.halfLength"));
    if (_config.getBool("stm.FOVcollimator.absorber.build",false)){
             finishNesting(collimatorFOVAbsorber,
                    findMaterialOrThrow(_config.getString("stm.FOVcollimator.absorber.material")),
                    0,
                    stmFOVCollAbsorberPositionInParent,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
    }

    if (verbosityLevel>0){
      std::cout<<__func__<<" STM FOV Coll (lead) z_center     = "<< stmFOVCollPositionInMu2e1.z() <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (lead) z_halflength = "<< stmFOVCollHalfLength1 <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (lead) z_min        = "<< stmFOVCollPositionInMu2e1.z()-stmFOVCollHalfLength1 <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (lead) z_max        = "<< stmFOVCollPositionInMu2e1.z()+stmFOVCollHalfLength1 <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) z_center     = "<< stmFOVCollPositionInMu2e2.z() <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) z_halflength = "<< stmFOVCollHalfLength2 <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) z_min        = "<< stmFOVCollPositionInMu2e2.z()-stmFOVCollHalfLength2 <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) z_max        = "<< stmFOVCollPositionInMu2e2.z()+stmFOVCollHalfLength2 <<std::endl;
      std::cout<<__func__<<" STM FOV Coll (poly) r_UpStr      = "<< pSTMFOVCollimatorParams.hole1RadiusUpStr() <<std::endl;
    }

    G4VSensitiveDetector *sd3 = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::STMFOVa());
    if(sd3) collimatorFOVAbsorber.logical->SetSensitiveDetector(sd3);

    //===================== Magnet and FOV Collimator Support Table ==========================

    //Just use a block of material for now (maybe stainless steel?, specified in configuration)
    G4Material*  stmMagnetSupportTableMaterial   = findMaterialOrThrow(pSTMMagnetSupportTableParams.materialName());
    const double stmMagnetSupportTableHalfLengths[3] = {pSTMMagnetSupportTableParams.tabletopHalfWidth(),
                                                 pSTMMagnetSupportTableParams.tabletopHalfHeight(),
                                                 pSTMMagnetSupportTableParams.tabletopHalfLength()};

    const double mstmMagnetStandLegRadius  = pSTMMagnetSupportTableParams.legRadius();

    //G4ThreeVector stmMagnetStandPositionInMu2e   = mstmMagnetPositionInMu2e   + G4ThreeVector(0.0, -(pSTMMagnetParams.zHalfLength()+stmMagnetSupportTableHalfLengths[1]), pipeDnStrExtentHalflength + 0.5*stmFOVCollUpStrSpace + stmFOVCollHalfLength1);
    //G4ThreeVector stmMagnetStandPositionInParent = mstmMagnetPositionInParent + G4ThreeVector(0.0, -(pSTMMagnetParams.zHalfLength()+stmMagnetSupportTableHalfLengths[1]), pipeDnStrExtentHalflength + 0.5*stmFOVCollUpStrSpace + stmFOVCollHalfLength1);
    G4ThreeVector stmMagnetSupportTablePositionInMu2e   = pSTMMagnetSupportTableParams.originInMu2e();
    G4ThreeVector stmMagnetSupportTablePositionInParent = pSTMMagnetSupportTableParams.originInMu2e() - parentCenterInMu2e;

    const double yExtentLow = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
    const double mstmMagnetStandLegHalfHeight = (yExtentLow-pSTMMagnetParams.zHalfLength()-2.0*stmMagnetSupportTableHalfLengths[1])/2.0;
    const TubsParams mstmMagnetStandLegParams(0.0, mstmMagnetStandLegRadius, mstmMagnetStandLegHalfHeight , 0.0, CLHEP::twopi);
    const double mstmMagnetStandLegOffsetX = stmMagnetSupportTableHalfLengths[0]  - mstmMagnetStandLegRadius - 1.0*CLHEP::cm;
    const double mstmMagnetStandLegOffsetZ = stmMagnetSupportTableHalfLengths[2] - mstmMagnetStandLegRadius - 1.0*CLHEP::cm;

    CLHEP::HepRotationX RXForLegs(90.0*CLHEP::degree);
    G4RotationMatrix *rotMatrixXforLegs = new G4RotationMatrix(RXForLegs);

    if (pSTMMagnetSupportTableParams.build()){
      VolumeInfo mstmMagnetStandInfo = nestBox("mstmMagnetStandPlatform",
                                          stmMagnetSupportTableHalfLengths,
                                          stmMagnetSupportTableMaterial,
                                          0x0,
                                          stmMagnetSupportTablePositionInParent, //mstmMagnetStandPositionInMother,
                                          parentInfo.logical,
                                          0,
                                          STMisVisible,
                                          G4Color::Gray(),
                                          STMisSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                         );

      VolumeInfo mstmMagnetStandLeg1Info = nestTubs( "mstmMagnetStandLeg1",
                                              mstmMagnetStandLegParams,
                                              stmMagnetSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmMagnetSupportTablePositionInParent+G4ThreeVector(mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight,-mstmMagnetStandLegOffsetZ),//mstmMagnetStandPositionInMother+G4ThreeVector(mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight,-mstmMagnetStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmMagnetStandLeg2Info = nestTubs( "mstmMagnetStandLeg2",
                                              mstmMagnetStandLegParams,
                                              stmMagnetSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmMagnetSupportTablePositionInParent+G4ThreeVector(mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight, mstmMagnetStandLegOffsetZ),//mstmMagnetStandPositionInMother+G4ThreeVector(mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight, mstmMagnetStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmMagnetStandLeg3Info = nestTubs( "mstmMagnetStandLeg3",
                                              mstmMagnetStandLegParams,
                                              stmMagnetSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmMagnetSupportTablePositionInParent+G4ThreeVector(-mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight,-mstmMagnetStandLegOffsetZ),//mstmMagnetStandPositionInMother+G4ThreeVector(-mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight,-mstmMagnetStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmMagnetStandLeg4Info = nestTubs( "mstmMagnetStandLeg4",
                                              mstmMagnetStandLegParams,
                                              stmMagnetSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmMagnetSupportTablePositionInParent+G4ThreeVector(-mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight, mstmMagnetStandLegOffsetZ),//mstmMagnetStandPositionInMother+G4ThreeVector(-mstmMagnetStandLegOffsetX,-stmMagnetSupportTableHalfLengths[1]-mstmMagnetStandLegHalfHeight, mstmMagnetStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

    }

		//======================Absorber to reduce event rate on STM=================
		const double stmAbsorberHalfLength1 = pSTMAbsorberParams.halfLength(); //wyq
		const double stmAbsorberHalfWidth1  = pSTMAbsorberParams.halfWidth();//wyq
		const double stmAbsorberHalfHeight1 = pSTMAbsorberParams.halfHeight();//wyq
		// position of absorber
		G4ThreeVector stmAbsorberPositionInMu2e1   = pSTMAbsorberParams.originInMu2e();//wyq
		G4ThreeVector stmAbsorberPositionInParent1 = pSTMAbsorberParams.originInMu2e() - parentCenterInMu2e;//wyq
		// Make the box for the absorber
		if (pSTMAbsorberParams.build()&&stmAbsorberHalfLength1>0.0){//wyq
			G4Box* boxAbsorber = new G4Box("boxAbsorber",stmAbsorberHalfWidth1,stmAbsorberHalfHeight1,stmAbsorberHalfLength1);//wyq
			VolumeInfo absorber;//wyq
			absorber.name = "absorber";//wyq
			absorber.solid=boxAbsorber; //wyq
			finishNesting(absorber,//wyq
					findMaterialOrThrow(pSTMAbsorberParams.material()),//wyq
					0,//wyq
					stmAbsorberPositionInParent1,//wyq
					parentInfo.logical,//wyq
					0,//wyq
					true,//wyq
					G4Colour::Magenta(),//wyq
					false,//wyq
					forceAuxEdgeVisible,//wyq
					placePV,//wyq
					doSurfaceCheck);//wyq
		}//wyq
		//======================Basorber as a scintilator STM=================
		const double stmBasorberHalfLength1 = pSTMBasorberParams.halfLength(); //wyq
		const double stmBasorberHalfWidth1  = pSTMBasorberParams.halfWidth();//wyq
		const double stmBasorberHalfHeight1 = pSTMBasorberParams.halfHeight();//wyq
		// position of stmFOVd1
		G4ThreeVector stmBasorberPositionInMu2e1   = pSTMBasorberParams.originInMu2e();//wyq
		G4ThreeVector stmBasorberPositionInParent1 = pSTMBasorberParams.originInMu2e() - parentCenterInMu2e;//wyq
		G4ThreeVector stmBasorberPositionInParent2 = stmBasorberPositionInParent1+G4ThreeVector(0,0,2*stmBasorberHalfLength1);//wyq
		// Make the box for the stmFOVd1
		if (pSTMBasorberParams.build()&&stmBasorberHalfLength1>0.0){//wyq
			G4Box* boxBasorber = new G4Box("boxBasorber",stmBasorberHalfWidth1,stmBasorberHalfHeight1,stmBasorberHalfLength1);//wyq
			G4Box* boxSTMFOVdd = new G4Box("boxSTMFOVdd",stmBasorberHalfWidth1,stmBasorberHalfHeight1,stmBasorberHalfLength1);//wyq
			VolumeInfo stmFOVd1;//wyq
			stmFOVd1.name = "stmFOVd1";//wyq
			stmFOVd1.solid=boxBasorber; //wyq
			VolumeInfo stmFOVd2;//wyq
			stmFOVd2.name = "stmFOVd2";//wyq
			stmFOVd2.solid=boxSTMFOVdd; //wyq
			finishNesting(stmFOVd1,//wyq
					findMaterialOrThrow(pSTMBasorberParams.material()),//wyq
					0,//wyq
					stmBasorberPositionInParent1,//wyq
					parentInfo.logical,//wyq
					0,//wyq
					true,//wyq
					G4Colour::Magenta(),//wyq
					false,//wyq
					forceAuxEdgeVisible,//wyq
					placePV,//wyq
					doSurfaceCheck);//wyq
      G4VSensitiveDetector *sd4 = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::STMFOVd());
      if(sd4) stmFOVd1.logical->SetSensitiveDetector(sd4);
			finishNesting(stmFOVd2,//wyq
					findMaterialOrThrow(pSTMBasorberParams.material()),//wyq
					0,//wyq
					stmBasorberPositionInParent2,//wyq
					parentInfo.logical,//wyq
					1,//wyq
					true,//wyq
					G4Colour::Magenta(),//wyq
					false,//wyq
					forceAuxEdgeVisible,//wyq
					placePV,//wyq
					doSurfaceCheck);//wyq
    G4VSensitiveDetector *sd5 = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::STMFOVd());
    if(sd5) stmFOVd2.logical->SetSensitiveDetector(sd5);
		}//wyq

    //===================== Spot-Size (SS) Collimator ==========================

    const double stmSSCollHalfLength1 = pSTMSSCollimatorParams.halfLength();
    const double stmSSCollHalfWidth1  = pSTMSSCollimatorParams.halfWidth();
    const double stmSSCollHalfHeight1 = pSTMSSCollimatorParams.halfHeight();
    const double stmSSCollHalfLength2 = pSTMSSCollimatorParams.linerHalfLength();
    const double stmSSCollHalfWidth2  = pSTMSSCollimatorParams.linerHalfWidth();
    const double stmSSCollHalfHeight2 = pSTMSSCollimatorParams.linerHalfHeight();
		const double stmSSCollLinerCutOutHalfLength=pSTMSSCollimatorParams.linerCutOutHalfLength();

    // position of SS collimator
    G4ThreeVector stmSSCollPositionInMu2e1   = pSTMSSCollimatorParams.originInMu2e();
    G4ThreeVector stmSSCollPositionInParent1 = pSTMSSCollimatorParams.originInMu2e() - parentCenterInMu2e+stmSSCollShift;
    // Make the box for the collimator
    G4Box* boxSSColl = new G4Box("boxSSColl",stmSSCollHalfWidth1,stmSSCollHalfHeight1,stmSSCollHalfLength1);
    //---
/*    //Make the tube(s) for the hole(s)
    G4Tubs *tubSSColl1 = new G4Tubs("tubSSColl1", 0.0, pSTMSSCollimatorParams.hole1RadiusUpStr(), stmSSCollHalfLength1+1.0, 0.0, CLHEP::twopi );
    G4Tubs *tubSSColl2 = new G4Tubs("tubSSColl2", 0.0, pSTMSSCollimatorParams.hole2RadiusUpStr(), stmSSCollHalfLength1+1.0, 0.0, CLHEP::twopi );
    G4SubtractionSolid *collimatorSStemp1 = new G4SubtractionSolid("collimatorSStemp1",boxSSColl,         tubSSColl1,0,G4ThreeVector(pSTMSSCollimatorParams.hole1xOffset(),0.0,0.0));
    G4SubtractionSolid *collimatorSStemp2 = 0;
    if (pSTMSSCollimatorParams.hole2Build()){
      collimatorSStemp2 = new G4SubtractionSolid("collimatorSStemp2",collimatorSStemp1, tubSSColl2,0,G4ThreeVector(pSTMSSCollimatorParams.hole2xOffset(),0.0,0.0));
    } else {
      collimatorSStemp2 = collimatorSStemp1;
    } */
    //---

    GeomHandle<StoppingTarget> stoppingTarget;
    TargetFoil const& foil_downstream = stoppingTarget->foil(stoppingTarget->nFoils()-1);
    const double z_tgtfoil_downstream = foil_downstream.centerInMu2e().z();
    const double z_collimator_downstream = stmSSCollPositionInMu2e1.z() + stmSSCollHalfLength1;
    const double z_distance_tgt_coll = (z_collimator_downstream-z_tgtfoil_downstream);

    //Make the conical.disk for the first hole as wide as the StoppingTarget+extra on one end
    //and as narrow as the desired collimation on the other end
		double collWindow1RmaxUp=foil_downstream.rOut()+150.0;//wyq
		double collWindow1RmaxDn=pSTMSSCollimatorParams.hole1RadiusDnStr();//wyq
    G4Tubs* _collWindow1=new G4Tubs("_collWindow1",//for a staight hole
				0.0,collWindow1RmaxDn,stmSSCollHalfLength2+0.01,0.0,CLHEP::twopi);
		//wyq if need liner absorber, enlarge the holes accordingly
		if(pSTMSSCollimatorParams.hole1Liner1Build()){//wyq
			collWindow1RmaxUp+=pSTMSSCollimatorParams.hole1Liner1Thickness();//wyq
			collWindow1RmaxDn+=pSTMSSCollimatorParams.hole1Liner1Thickness();//wyq
		}//wyq
		if(pSTMSSCollimatorParams.hole1Liner2Build()){//wyq
			collWindow1RmaxUp+=pSTMSSCollimatorParams.hole1Liner2Thickness();//wyq
			collWindow1RmaxDn+=pSTMSSCollimatorParams.hole1Liner2Thickness();//wyq
		}//wyq
    G4Cons* collWindow1 = new G4Cons( "collWindow1", //for a tilted hole
                                      0.0,                          // rMin cone upstream
				collWindow1RmaxUp, //wyq rMax cone upStream
				0.0,                          // rMin cone downstream
				collWindow1RmaxDn, //wyq rMax cone downstream
                                      z_distance_tgt_coll/2.0,      //halflength
                                      0.0,                          //start angle
                                      CLHEP::twopi                  //end angle
                                    );
    //Make the conical.disk for the second hole as wide as the StoppingTarget+extra on one end
    //and as narrow as the desired collimation on the other end
		double collWindow2RmaxUp=foil_downstream.rOut()+150.0;//wyq
		double collWindow2RmaxDn=pSTMSSCollimatorParams.hole2RadiusDnStr();//wyq
    G4Tubs* _collWindow2=new G4Tubs("_collWindow2",//for a staight hole
				0.0,collWindow2RmaxDn,stmSSCollHalfLength2+0.01,0.0,CLHEP::twopi);
		//wyq if need liner absorber, enlarge the holes accordingly
		if(pSTMSSCollimatorParams.hole2Liner1Build()){//wyq
			collWindow2RmaxUp+=pSTMSSCollimatorParams.hole2Liner1Thickness();//wyq
			collWindow2RmaxDn+=pSTMSSCollimatorParams.hole2Liner1Thickness();//wyq
		}//wyq
		if(pSTMSSCollimatorParams.hole2Liner2Build()){//wyq
			collWindow2RmaxUp+=pSTMSSCollimatorParams.hole2Liner2Thickness();//wyq
			collWindow2RmaxDn+=pSTMSSCollimatorParams.hole2Liner2Thickness();//wyq
		}//wyq
    G4Cons* collWindow2 = new G4Cons( "collWindow2",
                                      0.0,                          // rMin cone upstream
				collWindow2RmaxUp, //wyq rMax cone upStream
				0.0,                          // rMin cone downstream
				collWindow2RmaxDn, //wyq rMax cone downstream
                                      z_distance_tgt_coll/2.0,      //halflength
                                      0.0,                          //start angle
                                      CLHEP::twopi                  //end angle
                                    );


    const double xoffset_hole1 = pSTMSSCollimatorParams.hole1xOffset();
    const double angleY1 = -1.0*std::atan( (xoffset_hole1/2.0)/(z_distance_tgt_coll/2.0) );
    CLHEP::HepRotationY RYForCone1(angleY1);
    G4RotationMatrix *rotMatrixYforCone1 = new G4RotationMatrix(RYForCone1);
		//const double z_shift1 = z_distance_tgt_coll/2.0*std::sin(std::abs(angleY1)) + pSTMSSCollimatorParams.hole1RadiusDnStr()*std::sin(std::abs(angleY1));
		const double z_shift1 = z_distance_tgt_coll*2.0*std::sin(std::abs(angleY1/2.0))*std::sin(std::abs(angleY1/2.0)) + collWindow1RmaxDn*std::sin(std::abs(angleY1));//wyq make it precise

    const double xoffset_hole2 = pSTMSSCollimatorParams.hole2xOffset();
    const double angleY2 = -1.0*std::atan( (xoffset_hole2/2.0)/(z_distance_tgt_coll/2.0) );
    CLHEP::HepRotationY RYForCone2(angleY2);
    G4RotationMatrix *rotMatrixYforCone2 = new G4RotationMatrix(RYForCone2);
		//		const double z_shift2 = z_distance_tgt_coll/2.0*std::sin(std::abs(angleY2)) + pSTMSSCollimatorParams.hole2RadiusDnStr()*std::sin(std::abs(angleY2));
		const double z_shift2 = z_distance_tgt_coll*2.0*std::sin(std::abs(angleY2/2.0))*std::sin(std::abs(angleY2/2.0)) + collWindow2RmaxDn*std::sin(std::abs(angleY2));//wyq

//     const double z_FOVColl_downstream = stmFOVCollPositionInMu2e1.z() + stmFOVCollHalfLength1;
//     const double z_SScollimator_downstream = stmSSCollPositionInMu2e1.z() + stmSSCollHalfLength1;
//     const double z_distance_tgt_coll = (z_SScollimator_downstream-z_FOVColl_downstream);
//     G4Cons* collWindow1 = new G4Cons( "collWindow1",
//                                       0.0,                     // rMin upstream
//                                       pSTMFOVCollimatorParams.hole1RadiusUpStr()+5.0,  // rMax upStream, 0.5cm larger than FOV coll opening
//                                       0.0,                     // rMin downstream
//                                       pSTMSSCollimatorParams.hole1RadiusDnStr(), // rMax downstream
//                                       z_distance_tgt_coll/2.0, //halflength
//                                       0.0,                     //start angle
//                                       CLHEP::twopi             //end angle
//                                     );
//     //Make the conical.disk for the first hole as wide as the Stopping Target on one end
//     //and as narrow as the desired collimation on the other end
//     G4Cons* collWindow2 = new G4Cons( "collWindow2",
//                                       0.0,                     // rMin upstream
//                                       pSTMFOVCollimatorParams.hole1RadiusUpStr()+5.0,  // rMax upStream, 0.5cm larger than FOV coll opening
//                                       0.0,                     // rMin downstream
//                                       pSTMSSCollimatorParams.hole2RadiusUpStr(), // rMax downstream
//                                       z_distance_tgt_coll/2.0, //halflength
//                                       0.0,                     //start angle
//                                       CLHEP::twopi             //end angle
//                                     );

    // Combine into the Wall with the Hole
    // Use a G4SubtractionSolid to allow for another volume placement through it
    G4SubtractionSolid *collimatorSStemp1 = new G4SubtractionSolid("collimatorSStemp1",boxSSColl,collWindow1,rotMatrixYforCone1,G4ThreeVector(xoffset_hole1/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift1 ));
    G4SubtractionSolid *collimatorSStemp2 = 0;
    if (pSTMSSCollimatorParams.hole2Build()){
      collimatorSStemp2 = new G4SubtractionSolid("collimatorSStemp2",collimatorSStemp1,collWindow2,rotMatrixYforCone2,G4ThreeVector(xoffset_hole2/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift2 ));
    } else {
      collimatorSStemp2 = collimatorSStemp1;
    }

    //---

    //Make a box to subtract so liner can fit inside
    G4Box* boxSSCollLinerToSubt = new G4Box("boxSSCollLinerToSubt",stmSSCollHalfWidth2+2*stmSSCollLinerCutOutHalfLength+0.001,stmSSCollHalfHeight2+2*stmSSCollLinerCutOutHalfLength+0.001,stmSSCollHalfLength2+stmSSCollLinerCutOutHalfLength+0.001+100);
      G4Box* boxSSCollLinerToSubtLsub=new G4Box("boxSSCollLinerToSubtLsub",stmSSCollHalfWidth2/8,stmSSCollHalfHeight2+2*stmSSCollLinerCutOutHalfLength+0.001,stmSSCollHalfLength2*2/3);
      G4SubtractionSolid* boxSSCollLinerToSubtL=new G4SubtractionSolid("boxSSCollLinerToSubtL",boxSSCollLinerToSubt,boxSSCollLinerToSubtLsub,0,G4ThreeVector(-stmSSCollHalfWidth2*7/8-2*stmSSCollLinerCutOutHalfLength+0.001,0,-stmSSCollHalfLength2/3-stmSSCollLinerCutOutHalfLength-0.001));
    // Combine into the collimator with the liner cut-out and collimation hole
    VolumeInfo collimatorSS;
    collimatorSS.name = "collimatorSS";
    if (pSTMSSCollimatorParams.linerBuild()){
    if(stmSSCollShape==2){//L shape
       collimatorSS.solid = new G4SubtractionSolid(collimatorSS.name,collimatorSStemp2,boxSSCollLinerToSubtL,0,G4ThreeVector(0.0,0.0,0.0));
    }else{
       collimatorSS.solid = new G4SubtractionSolid(collimatorSS.name,collimatorSStemp2,boxSSCollLinerToSubt,0,G4ThreeVector(0.0,0.0,0.0));
    }
    } else {
       collimatorSS.solid = collimatorSStemp2;
    }

    //position of liner
    G4ThreeVector stmSSCollPositionInMu2e2   = stmSSCollPositionInMu2e1   + G4ThreeVector(0.0,0.0, -stmSSCollHalfLength1+stmSSCollHalfLength2);
    G4ThreeVector stmSSCollPositionInParent2 = stmSSCollPositionInParent1 + G4ThreeVector(0.0,0.0,0.0);// -stmSSCollHalfLength1+stmSSCollHalfLength2);
    // make the box for the liner
    G4Box* boxSSCollLiner = 0;
    G4SubtractionSolid *collimatorSSLinerTemp1 = 0;
    if (pSTMSSCollimatorParams.linerBuild()){
      boxSSCollLiner = new G4Box("boxSSCollLiner",stmSSCollHalfWidth2,stmSSCollHalfHeight2,stmSSCollHalfLength2);
      if(stmSSCollxshift==0&&stmSSCollyshift==0)collimatorSSLinerTemp1 = new G4SubtractionSolid("collimatorSSLinerTemp1",boxSSCollLiner,collWindow1,rotMatrixYforCone1,G4ThreeVector(xoffset_hole1/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift1 ));
			else collimatorSSLinerTemp1 = new G4SubtractionSolid("collimatorSSLinerTemp1",boxSSCollLiner,_collWindow1,0,G4ThreeVector(xoffset_hole1,0.0,0.0)-stmSSCollShift);
    } else {
      collimatorSSLinerTemp1 = 0;
    }
    // tubSSColl1,0,(stmSSCollPositionInMu2e1-stmSSCollPositionInMu2e1)+G4ThreeVector(pSTMSSCollimatorParams.hole1xOffset(),0.0,0.0));
    G4SubtractionSolid *collimatorSSLinerTemp2 = 0;
    if (pSTMSSCollimatorParams.hole2Build()){
      if(stmSSCollxshift==0&&stmSSCollyshift==0)collimatorSSLinerTemp2 = new G4SubtractionSolid("collimatorSSLinerTemp2",collimatorSSLinerTemp1,collWindow2,rotMatrixYforCone2,G4ThreeVector(xoffset_hole2/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift2 ));
			else collimatorSSLinerTemp2 = new G4SubtractionSolid("collimatorSSLinerTemp2",collimatorSSLinerTemp1,_collWindow2,0,G4ThreeVector(xoffset_hole2,0.0,0.0)-stmSSCollShift);
      //tubSSColl2,0,(stmSSCollPositionInMu2e1-stmSSCollPositionInMu2e1)+G4ThreeVector(pSTMSSCollimatorParams.hole2xOffset(),0.0,0.0));
    } else {
      collimatorSSLinerTemp2 = collimatorSSLinerTemp1;
    }
    G4Box *collimatorSSLinerTemp2Lsub=new G4Box("collimatorSSLinerTemp2Lsub",stmSSCollHalfWidth2/8,stmSSCollHalfHeight2+0.001,stmSSCollHalfLength2*2/3);
    G4SubtractionSolid *collimatorSSLinerTemp2L=new G4SubtractionSolid("collimatorSSLinerTemp2L",collimatorSSLinerTemp2,collimatorSSLinerTemp2Lsub,0,G4ThreeVector(-stmSSCollHalfWidth2*7/8-0.001,0,-stmSSCollHalfLength2/3-0.001));

    VolumeInfo collimatorSSliner;
    collimatorSSliner.name = "collimatorSSliner";
    if(stmSSCollShape==2){
      collimatorSSliner.solid = collimatorSSLinerTemp2L;
    }else{
      collimatorSSliner.solid = collimatorSSLinerTemp2;
    }
		/////__________________//wyq
		//wyq hole1 liner1
		G4Cons* hole1Liner1 = new G4Cons( "hole1Liner1", //wyq
				collWindow1RmaxUp-pSTMSSCollimatorParams.hole1Liner1Thickness(),                          // rMin cone upstream//wyq
				collWindow1RmaxUp, //wyq rMax cone upStream//wyq
				collWindow1RmaxDn-pSTMSSCollimatorParams.hole1Liner1Thickness(),                          // rMin cone downstream//wyq
				collWindow1RmaxDn, //wyq rMax cone nstream//wyq
				z_distance_tgt_coll/2.0,      //halflength//wyq
				0.0,                          //start //wyq
				CLHEP::twopi                  //end angle //wyq
				);	//wyq
		G4IntersectionSolid *hole1Liner1a=new G4IntersectionSolid("hole1Liner1a",boxSSColl,hole1Liner1,rotMatrixYforCone1,G4ThreeVector(xoffset_hole1/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift1 ));  //wyq
		//wyq hole1 liner 2
		G4Cons* hole1Liner2 = new G4Cons( "hole1Liner2", //wyq
				collWindow1RmaxUp-pSTMSSCollimatorParams.hole1Liner1Thickness()-pSTMSSCollimatorParams.hole1Liner2Thickness(),                          // rMin cone upstream//wyq
				collWindow1RmaxUp-pSTMSSCollimatorParams.hole1Liner1Thickness(), //wyq rMax cone upStream//wyq
				collWindow1RmaxDn-pSTMSSCollimatorParams.hole1Liner1Thickness()-pSTMSSCollimatorParams.hole1Liner2Thickness(),                          // rMin cone downstream//wyq
				collWindow1RmaxDn-pSTMSSCollimatorParams.hole1Liner1Thickness(), //wyq rMax cone nstream//wyq
				z_distance_tgt_coll/2.0,      //halflength//wyq
				0.0,                          //start //wyq
				CLHEP::twopi                  //end angle //wyq
				);	//wyq
		G4IntersectionSolid *hole1Liner2a=new G4IntersectionSolid("hole1Liner2a",boxSSColl,hole1Liner2,rotMatrixYforCone1,G4ThreeVector(xoffset_hole1/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift1 ));  //wyq
		//wyq hole2 liner1
		G4Cons* hole2Liner1 = new G4Cons( "hole2Liner1", //wyq
				collWindow2RmaxUp-pSTMSSCollimatorParams.hole2Liner1Thickness(),                          // rMin cone upstream//wyq
				collWindow2RmaxUp, //wyq rMax cone upStream//wyq
				collWindow2RmaxDn-pSTMSSCollimatorParams.hole2Liner1Thickness(),                          // rMin cone downstream//wyq
				collWindow2RmaxDn, //wyq rMax cone nstream//wyq
				z_distance_tgt_coll/2.0,      //halflength//wyq
				0.0,                          //start //wyq
				CLHEP::twopi                  //end angle //wyq
				);	//wyq
		G4IntersectionSolid *hole2Liner1a=new G4IntersectionSolid("hole2Liner1a",boxSSColl,hole2Liner1,rotMatrixYforCone2,G4ThreeVector(xoffset_hole2/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift2 ));  //wyq
		//wyq hole2 liner 2
		G4Cons* hole2Liner2 = new G4Cons( "hole2Liner2", //wyq
				collWindow2RmaxUp-pSTMSSCollimatorParams.hole2Liner1Thickness()-pSTMSSCollimatorParams.hole2Liner2Thickness(),                          // rMin cone upstream//wyq
				collWindow2RmaxUp-pSTMSSCollimatorParams.hole2Liner1Thickness(), //wyq rMax cone upStream//wyq
				collWindow2RmaxDn-pSTMSSCollimatorParams.hole2Liner1Thickness()-pSTMSSCollimatorParams.hole2Liner2Thickness(),                          // rMin cone downstream//wyq
				collWindow2RmaxDn-pSTMSSCollimatorParams.hole2Liner1Thickness(), //wyq rMax cone nstream//wyq
				z_distance_tgt_coll/2.0,      //halflength//wyq
				0.0,                          //start //wyq
				CLHEP::twopi                  //end angle //wyq
				);	//wyq
		G4IntersectionSolid *hole2Liner2a=new G4IntersectionSolid("hole2Liner2a",boxSSColl,hole2Liner2,rotMatrixYforCone2,G4ThreeVector(xoffset_hole2/2.0,0.0,-1.0*z_distance_tgt_coll/2.0+stmSSCollHalfLength1+z_shift2 ));  //wyq

		VolumeInfo hole1Liner1b;//wyq
		hole1Liner1b.name = "hole1Liner1b";//wyq
		hole1Liner1b.solid=hole1Liner1a; //wyq
		VolumeInfo hole1Liner2b;//wyq
		hole1Liner2b.name = "hole1Liner2b";//wyq
		hole1Liner2b.solid=hole1Liner2a; //wyq

		VolumeInfo hole2Liner1b;//wyq
		hole2Liner1b.name = "hole2Liner1b";//wyq
		hole2Liner1b.solid=hole2Liner1a; //wyq
		VolumeInfo hole2Liner2b;//wyq
		hole2Liner2b.name = "hole2Liner2b";//wyq
		hole2Liner2b.solid=hole2Liner2a; //wyq
		///_______________wyq
    if (pSTMSSCollimatorParams.build()){
      finishNesting(collimatorSS,
                    findMaterialOrThrow(pSTMSSCollimatorParams.material()),
                    0,
                    stmSSCollPositionInParent1,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
      if (pSTMSSCollimatorParams.linerBuild()){
        finishNesting(collimatorSSliner,
                      findMaterialOrThrow(pSTMSSCollimatorParams.linerMaterial()),
                      0,
                      stmSSCollPositionInParent2,
                      parentInfo.logical,
                      0,
                      STMisVisible,
                      G4Colour::Magenta(),
                      STMisSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck);
      }
			if(pSTMSSCollimatorParams.hole1Liner1Build()){//wyq
				finishNesting(hole1Liner1b,//wyq
						findMaterialOrThrow(pSTMSSCollimatorParams.hole1Liner1Material()),//wyq
						rotMatrixYforCone1,//wyq
						stmSSCollPositionInParent2,//wyq
						parentInfo.logical,//wyq
						0,true,G4Color::Green(),false,forceAuxEdgeVisible,placePV,doSurfaceCheck);//wyq
			}//wyq
			if(pSTMSSCollimatorParams.hole1Liner2Build()){//wyq
				finishNesting(hole1Liner2b,//wyq
						findMaterialOrThrow(pSTMSSCollimatorParams.hole1Liner2Material()),//wyq
						rotMatrixYforCone1,//wyq
						stmSSCollPositionInParent2,//wyq
						parentInfo.logical,//wyq
						0,true,G4Color::Red(),false,forceAuxEdgeVisible,placePV,doSurfaceCheck);//wyq
			}//wyq
			if(pSTMSSCollimatorParams.hole2Liner1Build()){//wyq
				finishNesting(hole2Liner1b,//wyq
						findMaterialOrThrow(pSTMSSCollimatorParams.hole2Liner1Material()),//wyq
						rotMatrixYforCone1,//wyq
						stmSSCollPositionInParent2,//wyq
						parentInfo.logical,//wyq
						0,true,G4Colour::Green(),false,forceAuxEdgeVisible,placePV,doSurfaceCheck);//wyq
			}//wyq
			if(pSTMSSCollimatorParams.hole2Liner2Build()){//wyq
				finishNesting(hole2Liner2b,//wyq
						findMaterialOrThrow(pSTMSSCollimatorParams.hole2Liner2Material()),//wyq
						rotMatrixYforCone1,//wyq
						stmSSCollPositionInParent2,//wyq
						parentInfo.logical,//wyq
						0,true,G4Colour::Green(),false,forceAuxEdgeVisible,placePV,doSurfaceCheck);//wyq
    }
		}

    if (verbosityLevel>0){
      std::cout<<__func__<<" STM SS Coll (lead)   z_center     = "<< stmSSCollPositionInMu2e1.z() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)   z_halflength = "<< stmSSCollHalfLength1 <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)   z_max        = "<< stmSSCollPositionInMu2e1.z()+stmSSCollHalfLength1 <<std::endl;
      std::cout<<__func__<<" STM SS Coll (lead)   r_DnStr      = "<< pSTMSSCollimatorParams.hole1RadiusDnStr() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (cutout) z_center     = "<< stmSSCollPositionInMu2e2.z() <<std::endl;
      std::cout<<__func__<<" STM SS Coll (cutout) z_halflength = "<< stmSSCollHalfLength2 <<std::endl;
      std::cout<<__func__<<" STM SS Coll (cutout) z_min        = "<< stmSSCollPositionInMu2e2.z()-stmSSCollHalfLength2 <<std::endl;
    }



    //===================== STM Detector Support Table ==========================

    //Just use a block of material for now (maybe stainless steel?, specified in configuration)
    G4Material*  stmDetectorSupportTableMaterial   = findMaterialOrThrow(pSTMDetectorSupportTableParams.materialName());
    const double stmDetectorSupportTableHalfLengths[3] = {pSTMDetectorSupportTableParams.tabletopHalfWidth(),
                                                 pSTMDetectorSupportTableParams.tabletopHalfHeight(),
                                                 pSTMDetectorSupportTableParams.tabletopHalfLength()};

    const double mstmDetectorStandLegRadius  = pSTMDetectorSupportTableParams.legRadius();

    G4ThreeVector stmDetectorSupportTablePositionInMu2e   = pSTMDetectorSupportTableParams.originInMu2e();
    G4ThreeVector stmDetectorSupportTablePositionInParent = pSTMDetectorSupportTableParams.originInMu2e() - parentCenterInMu2e+G4ThreeVector(0,-STMShield_thick_bot,0.0);

    //const double yExtentLow = std::abs(_config.getDouble("yOfFloorSurface.below.mu2eOrigin") );
    const double mstmDetectorStandLegHalfHeight = (yExtentLow-pSTMSSCollimatorParams.halfHeight()-2.0*stmDetectorSupportTableHalfLengths[1])/2.0-STMShield_thick_bot/2.0;
    const TubsParams mstmDetectorStandLegParams(0.0, mstmDetectorStandLegRadius, mstmDetectorStandLegHalfHeight , 0.0, CLHEP::twopi);
    const double mstmDetectorStandLegOffsetX = stmDetectorSupportTableHalfLengths[0]  - mstmDetectorStandLegRadius - 1.0*CLHEP::cm;
    const double mstmDetectorStandLegOffsetZ = stmDetectorSupportTableHalfLengths[2] - mstmDetectorStandLegRadius - 1.0*CLHEP::cm;

    //CLHEP::HepRotationX RXForLegs(90.0*CLHEP::degree);
    //G4RotationMatrix *rotMatrixXforLegs = new G4RotationMatrix(RXForLegs);

    if (pSTMDetectorSupportTableParams.build()){
      VolumeInfo mstmDetectorStandInfo = nestBox("mstmDetectorStandPlatform",
                                          stmDetectorSupportTableHalfLengths,
                                          stmDetectorSupportTableMaterial,
                                          0x0,
                                          stmDetectorSupportTablePositionInParent, //mstmDetectorStandPositionInMother,
                                          parentInfo.logical,
                                          0,
                                          STMisVisible,
                                          G4Color::Gray(),
                                          STMisSolid,
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                         );

      VolumeInfo mstmDetectorStandLeg1Info = nestTubs( "mstmDetectorStandLeg1",
                                              mstmDetectorStandLegParams,
                                              stmDetectorSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmDetectorSupportTablePositionInParent+G4ThreeVector(mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight,-mstmDetectorStandLegOffsetZ),//mstmDetectorStandPositionInMother+G4ThreeVector(mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight,-mstmDetectorStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmDetectorStandLeg2Info = nestTubs( "mstmDetectorStandLeg2",
                                              mstmDetectorStandLegParams,
                                              stmDetectorSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmDetectorSupportTablePositionInParent+G4ThreeVector(mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight, mstmDetectorStandLegOffsetZ),//mstmDetectorStandPositionInMother+G4ThreeVector(mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight, mstmDetectorStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmDetectorStandLeg3Info = nestTubs( "mstmDetectorStandLeg3",
                                              mstmDetectorStandLegParams,
                                              stmDetectorSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmDetectorSupportTablePositionInParent+G4ThreeVector(-mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight,-mstmDetectorStandLegOffsetZ),//mstmDetectorStandPositionInMother+G4ThreeVector(-mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight,-mstmDetectorStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

      VolumeInfo mstmDetectorStandLeg4Info = nestTubs( "mstmDetectorStandLeg4",
                                              mstmDetectorStandLegParams,
                                              stmDetectorSupportTableMaterial,
                                              rotMatrixXforLegs,
                                              stmDetectorSupportTablePositionInParent+G4ThreeVector(-mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight, mstmDetectorStandLegOffsetZ),//mstmDetectorStandPositionInMother+G4ThreeVector(-mstmDetectorStandLegOffsetX,-stmDetectorSupportTableHalfLengths[1]-mstmDetectorStandLegHalfHeight, mstmDetectorStandLegOffsetZ),
                                              parentInfo,
                                              0,
                                              STMisVisible,
                                              G4Color::Gray(),
                                              STMisSolid,
                                              forceAuxEdgeVisible,
                                              placePV,
                                              doSurfaceCheck
                                              );

    }


    //===================== STM Detector 1 ==========================
///   if(pSTMDetector1Params.build()){
///    G4Material*  stmDet1Material                 =  findMaterialOrThrow(pSTMDetector1Params.crystalMaterial());
///    const double stmDet1ROut                     =  pSTMDetector1Params.crystalRadiusOut();
///    const double stmDet1HalfLength               =  pSTMDetector1Params.crystalHalfLength();
///    G4Material*  stmDet1CanMaterial              =  findMaterialOrThrow(pSTMDetector1Params.canMaterial());
///    const double stmDet1CanRIn                   =  pSTMDetector1Params.canRadiusIn();
///    const double stmDet1CanROut                  =  pSTMDetector1Params.canRadiusOut();
///    const double stmDet1CanHalfLength            =  pSTMDetector1Params.canHalfLength();
///    G4Material*  stmDet1CanUpStrWindowMaterial   =  findMaterialOrThrow(pSTMDetector1Params.canUpStrWindowMaterial());
///    const double stmDet1CanUpStrWindowHalfLength =  pSTMDetector1Params.canUpStrWindowHalfLength();
///    G4Material*  stmDet1CanDnStrWindowMaterial   =  stmDet1CanMaterial;
///    const double stmDet1CanDnStrWindowHalfLength =  stmDet1CanROut - stmDet1CanRIn;
///
///    if ( stmDet1ROut > stmDet1CanRIn ){
///      throw cet::exception("GEOM")<< " STM: det1 radius is larger than the inner radius of the can. \n" ;
///    }
///    if ( stmDet1HalfLength > stmDet1CanHalfLength-stmDet1CanDnStrWindowHalfLength-stmDet1CanUpStrWindowHalfLength ){
///      throw cet::exception("GEOM")<< " STM: det1 crystal length is larger than the inner length of the can. \n" ;
///    }
///
///    const TubsParams stmDet1Params(0., stmDet1ROut, stmDet1HalfLength);
///    const TubsParams stmDet1CanParams(stmDet1CanRIn, stmDet1CanROut, stmDet1CanHalfLength);
///    const TubsParams stmDet1CanGasParams(0.,  stmDet1CanRIn, stmDet1CanHalfLength - stmDet1CanUpStrWindowHalfLength - stmDet1CanDnStrWindowHalfLength);
///    const TubsParams stmDet1CanUpStrWindowParams(0., stmDet1CanRIn, stmDet1CanUpStrWindowHalfLength);
///    const TubsParams stmDet1CanDnStrWindowParams(0., stmDet1CanRIn, stmDet1CanDnStrWindowHalfLength);
///
///    G4ThreeVector stmDet1CanPositionInMu2e   = pSTMDetector1Params.originInMu2e();
///    G4ThreeVector stmDet1CanPositionInParent = pSTMDetector1Params.originInMu2e() - parentCenterInMu2e;
///    G4ThreeVector stmDet1PositionInMu2e      = stmDet1CanPositionInMu2e;
///    G4ThreeVector stmDet1PositionInParent    = stmDet1CanPositionInParent;
///
///    VolumeInfo stmDet1CanInfo = nestTubs( "stmDet1Can",
///                                         stmDet1CanParams,
///                                         stmDet1CanMaterial,
///                                         0x0,
///                                         stmDet1CanPositionInParent,
///                                         parentInfo,
///                                         0,
///                                         STMisVisible,
///                                         G4Color::Green(),
///                                         STMisSolid,
///                                         forceAuxEdgeVisible,
///                                         placePV,
///                                         doSurfaceCheck
///                                         );
///
///    VolumeInfo stmDet1 = nestTubs("stmDet1",
///                                      stmDet1Params,
///                                      stmDet1Material,
///                                      0x0,
///                                      stmDet1CanPositionInParent,
///                                      parentInfo,
///                                      0,
///                                      STMisVisible,
///                                      G4Color::Red(),
///                                      STMisSolid,
///                                      forceAuxEdgeVisible,
///                                      placePV,
///                                      doSurfaceCheck
///                                      );
///
///    VolumeInfo stmDet1CanUpStrWindowInfo = nestTubs( "stmDet1CanUpStrWindow",
///                                            stmDet1CanUpStrWindowParams,
///                                            stmDet1CanUpStrWindowMaterial,
///                                            0x0,
///                                            stmDet1CanPositionInParent + G4ThreeVector(0.0,0.0,-1.0*stmDet1CanHalfLength + stmDet1CanUpStrWindowHalfLength),
///                                            parentInfo,
///                                            0,
///                                            STMisVisible,
///                                            G4Color::Green(),
///                                            STMisSolid,
///                                            forceAuxEdgeVisible,
///                                            placePV,
///                                            doSurfaceCheck
///                                            );
///
///    VolumeInfo stmDet1CanDnStrWindowInfo = nestTubs( "stmDet1CanDnStrWindow",
///                                            stmDet1CanDnStrWindowParams,
///                                            stmDet1CanDnStrWindowMaterial,
///                                            0x0,
///                                            stmDet1CanPositionInParent + G4ThreeVector(0.0,0.0, stmDet1CanHalfLength - stmDet1CanDnStrWindowHalfLength),
///                                            parentInfo,
///                                            0,
///                                            STMisVisible,
///                                            G4Color::Green(),
///                                            STMisSolid,
///                                            forceAuxEdgeVisible,
///                                            placePV,
///                                            doSurfaceCheck
///                                            );
///
///    if (verbosityLevel>0){
///      std::cout << __func__ << " Warning: Gas not implemented inside STM detector1 can! (so that VD inside can does not overlap with can gas)" << std::endl;
///    }
///    // Make stmDet1 a sensitive detector.
///    G4VSensitiveDetector *sd1 = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::STMDet());
///    if(sd1) stmDet1.logical->SetSensitiveDetector(sd1);
///
///	 }

    //===================== STM Detector 2 ==========================
//>===shielding
//double STMShield_thick_front=2*STMShield_Ttop_Poly-STMShield_Ttop_Pb2-STMShield_Ttop_Cu-STMShield_Tbot_Al;
//double STMShield_Xin1=STMShield_X1-STMShield_thick-STMShield_thick/fabs(sin(STMShield_right_Theta));
//double STMShield_Xin2=STMShield_X2-STMShield_thick-STMShield_thick/fabs(sin(STMShield_right_Theta))-STMShield_Zgap*tan(fabs(STMShield_right_Theta));
//double STM_Shield_Zshift=110;
G4Trap* STMShield_top_Al=new G4Trap("STMShield_top_Al",0.5*STMShield_Z1,STMShield_Theta,0,0.5*STMShield_Ttop_Al,0.5*STMShield_X1,0.5*STMShield_X1,0,0.5*STMShield_Ttop_Al,0.5*STMShield_X2,0.5*STMShield_X2,0);
G4Trap* STMShield_bot_Al_1=new G4Trap("STMShield_bot_Al_1",0.5*STMShield_Z1,STMShield_Theta,0,0.5*STMShield_Tbot_Al,0.5*STMShield_X1,0.5*STMShield_X1,0,0.5*STMShield_Tbot_Al,0.5*STMShield_X2,0.5*STMShield_X2,0);
G4Box* STMShield_bot_Al_2=new G4Box("STMShield_bot_Al_2",stmSSCollHalfWidth1,0.5*STMShield_Tbot_Al,stmSSCollHalfLength1);
G4UnionSolid* STMShield_bot_Al=new G4UnionSolid("STMShield_bot_Al",STMShield_bot_Al_1,STMShield_bot_Al_2,0,G4ThreeVector(STMShield_Ttop_CenterX,0.0,-0.5*STMShield_Z1-0.5*STMShield_thick_front));
G4Trap* STMShield_top_Cu=new G4Trap("STMShield_top_Cu",0.5*STMShield_Z1,STMShield_Theta,0,0.5*STMShield_Ttop_Cu,0.5*STMShield_X1,0.5*STMShield_X1,0,0.5*STMShield_Ttop_Cu,0.5*STMShield_X2,0.5*STMShield_X2,0);
G4Box* STMShield_bot_Cu_2=new G4Box("STMShield_bot_Cu_2",stmSSCollHalfWidth1,0.5*STMShield_Ttop_Cu,stmSSCollHalfLength1);
G4UnionSolid* STMShield_bot_Cu=new G4UnionSolid("STMShield_bot_Cu",STMShield_top_Cu,STMShield_bot_Cu_2,0,G4ThreeVector(STMShield_Ttop_CenterX,0.0,-0.5*STMShield_Z1-0.5*STMShield_thick_front));
G4Trap* STMShield_top_Pb1=new G4Trap("STMShield_top_Pb1",0.5*STMShield_Z1,STMShield_Theta,0,0.5*STMShield_Ttop_Pb1,0.5*STMShield_X1,0.5*STMShield_X1,0,0.5*STMShield_Ttop_Pb1,0.5*STMShield_X2,0.5*STMShield_X2,0);
G4Trap* STMShield_top_Pb2=new G4Trap("STMShield_top_Pb2",0.5*STMShield_Z1,STMShield_Theta,0,0.5*STMShield_Ttop_Pb2,0.5*STMShield_X1,0.5*STMShield_X1,0,0.5*STMShield_Ttop_Pb2,0.5*STMShield_X2,0.5*STMShield_X2,0);
G4Box* STMShield_bot_Pb2_2=new G4Box("STMShield_bot_Pb2_2",stmSSCollHalfWidth1,0.5*STMShield_Ttop_Pb2,stmSSCollHalfLength1);
G4UnionSolid* STMShield_bot_Pb2=new G4UnionSolid("STMShield_bot_Pb2",STMShield_top_Pb2,STMShield_bot_Pb2_2,0,G4ThreeVector(STMShield_Ttop_CenterX,0.0,-0.5*STMShield_Z1-0.5*STMShield_thick_front));
G4Trap* STMShield_top_Poly=new G4Trap("STMShield_top_Poly",0.5*STMShield_Z1,STMShield_Theta,0,0.5*STMShield_Ttop_Poly,0.5*STMShield_X1,0.5*STMShield_X1,0,0.5*STMShield_Ttop_Poly,0.5*STMShield_X2,0.5*STMShield_X2,0);
G4Box* STMShield_bot_Poly_2=new G4Box("STMShield_bot_Poly_2",stmSSCollHalfWidth1,0.5*STMShield_Ttop_Poly,stmSSCollHalfLength1);
G4UnionSolid* STMShield_bot_Poly=new G4UnionSolid("STMShield_bot_Poly",STMShield_top_Poly,STMShield_bot_Poly_2,0,G4ThreeVector(STMShield_Ttop_CenterX,0.0,-0.5*STMShield_Z1-0.5*STMShield_thick_front));

			G4ThreeVector stmShieldingPositionInMu2e=pSTMShieldingParams.originInMu2e();//G4ThreeVector(-3904,0,40900-stmShieldingHalfLength);//wyq
			G4ThreeVector stmShieldingPositionInParent=stmShieldingPositionInMu2e-parentCenterInMu2e;
			G4ThreeVector stmShieldingPositionInParentShiftX_top=stmShieldingPositionInParent+G4ThreeVector(-STMShield_Ttop_CenterX,STMShield_Y1/2.0+0.5*STMShield_Ttop_Al+STMShield_Ttop_Ycrack,0.0)+stmSSCollShift;//wyq
			G4ThreeVector stmShieldingPositionInParentShiftX_bot=stmShieldingPositionInParent+G4ThreeVector(-STMShield_Ttop_CenterX,-STMShield_Y1/2.0-0.5*STMShield_Tbot_Al,0.0)+stmSSCollShift;//wyq
    VolumeInfo stmShielding_top_1;
		stmShielding_top_1.name="stmShielding_top_1";
		stmShielding_top_1.solid=STMShield_top_Al;
		finishNesting(stmShielding_top_1,findMaterialOrThrow("G4_Al"),0,stmShieldingPositionInParentShiftX_top,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_top_2;
		stmShielding_top_2.name="stmShielding_top_2";
		stmShielding_top_2.solid=STMShield_top_Cu;
		finishNesting(stmShielding_top_2,findMaterialOrThrow("G4_Cu"),0,stmShieldingPositionInParentShiftX_top+G4ThreeVector(0,(STMShield_Ttop_Al+STMShield_Ttop_Cu+0.001)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_top_3;
		stmShielding_top_3.name="stmShielding_top_3";
		stmShielding_top_3.solid=STMShield_top_Pb2;
		finishNesting(stmShielding_top_3,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_top+G4ThreeVector(0,(STMShield_Ttop_Al+2*STMShield_Ttop_Cu+STMShield_Ttop_Pb2+0.002)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_top_4;
		stmShielding_top_4.name="stmShielding_top_4";
		stmShielding_top_4.solid=STMShield_top_Poly;
		finishNesting(stmShielding_top_4,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_top+G4ThreeVector(0,(STMShield_Ttop_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+STMShield_Ttop_Poly+0.003)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_top_5;
		stmShielding_top_5.name="stmShielding_top_5";
		stmShielding_top_5.solid=STMShield_top_Pb1;
		finishNesting(stmShielding_top_5,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_top+G4ThreeVector(0,(STMShield_Ttop_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly+STMShield_Ttop_Pb1+0.004)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_top_6;
		stmShielding_top_6.name="stmShielding_top_6";
		stmShielding_top_6.solid=STMShield_top_Poly;
		finishNesting(stmShielding_top_6,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_top+G4ThreeVector(0,(STMShield_Ttop_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly+2*STMShield_Ttop_Pb1+STMShield_Ttop_Poly+0.005)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
		//for bottom
    VolumeInfo stmShielding_bot_1;
		stmShielding_bot_1.name="stmShielding_bot_1";
		stmShielding_bot_1.solid=STMShield_bot_Al;
		finishNesting(stmShielding_bot_1,findMaterialOrThrow("G4_Al"),0,stmShieldingPositionInParentShiftX_bot,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_bot_2;
		stmShielding_bot_2.name="stmShielding_bot_2";
		stmShielding_bot_2.solid=STMShield_bot_Cu;
		finishNesting(stmShielding_bot_2,findMaterialOrThrow("G4_Cu"),0,stmShieldingPositionInParentShiftX_bot+G4ThreeVector(0,-(STMShield_Tbot_Al+STMShield_Ttop_Cu+0.001)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_bot_3;
		stmShielding_bot_3.name="stmShielding_bot_3";
		stmShielding_bot_3.solid=STMShield_bot_Pb2;
		finishNesting(stmShielding_bot_3,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_bot+G4ThreeVector(0,-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+STMShield_Ttop_Pb2+0.002)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_bot_4;
		stmShielding_bot_4.name="stmShielding_bot_4";
		stmShielding_bot_4.solid=STMShield_bot_Poly;
		finishNesting(stmShielding_bot_4,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_bot+G4ThreeVector(0,-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+STMShield_Ttop_Poly+0.003)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_bot_5;
		stmShielding_bot_5.name="stmShielding_bot_5";
		stmShielding_bot_5.solid=STMShield_bot_Pb2;
		finishNesting(stmShielding_bot_5,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_bot+G4ThreeVector(0,-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+0.004)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_bot_6;
		stmShielding_bot_6.name="stmShielding_bot_6";
		stmShielding_bot_6.solid=STMShield_bot_Poly;
		finishNesting(stmShielding_bot_6,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_bot+G4ThreeVector(0,-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly+2*STMShield_Ttop_Pb2+STMShield_Ttop_Poly+0.005)/2.0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
//for left
G4Box *STMShield_left_Al=new G4Box("STMShield_left_Al",0.5*STMShield_Tbot_Al,0.5*STMShield_Y1-0.001,0.5*STMShield_Z2);
G4Box *STMShield_left_Cu=new G4Box("STMShield_left_Cu",0.5*STMShield_Ttop_Cu,0.5*STMShield_Y1-0.001,0.5*STMShield_Z2);
G4Box *STMShield_left_Pb2=new G4Box("STMShield_left_Pb2",0.5*STMShield_Ttop_Pb2,0.5*STMShield_Y1-0.001,0.5*STMShield_Z2);
G4Box *STMShield_left_Poly=new G4Box("STMShield_left_Poly",0.5*STMShield_Ttop_Poly,0.5*STMShield_Y1-0.001,0.5*STMShield_Z2);
double STMShield_Tleft_CenterX=STMShield_X1/2.0-0.5*STMShield_Ttop_Al-STMShield_Ttop_Cu-STMShield_Ttop_Pb2*2-STMShield_Ttop_Poly*2;
G4ThreeVector stmShieldingPositionInParentShiftX_left=stmShieldingPositionInParent+G4ThreeVector(STMShield_Tleft_CenterX,0.0,0.5*(STMShield_Tbot_Al+2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Cu)-0.5*STMShield_Zgap)+stmSSCollShift;//wyq
    VolumeInfo stmShielding_left_1;
		stmShielding_left_1.name="stmShielding_left_1";
		stmShielding_left_1.solid=STMShield_left_Al;
		finishNesting(stmShielding_left_1,findMaterialOrThrow("G4_Al"),0,stmShieldingPositionInParentShiftX_left,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_left_2;
		stmShielding_left_2.name="stmShielding_left_2";
		stmShielding_left_2.solid=STMShield_left_Cu;
		finishNesting(stmShielding_left_2,findMaterialOrThrow("G4_Cu"),0,stmShieldingPositionInParentShiftX_left+G4ThreeVector((STMShield_Tbot_Al+STMShield_Ttop_Cu+0.001)/2.0,0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_left_3;
		stmShielding_left_3.name="stmShielding_left_3";
		stmShielding_left_3.solid=STMShield_left_Pb2;
		finishNesting(stmShielding_left_3,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_left+G4ThreeVector((STMShield_Tbot_Al+2*STMShield_Ttop_Cu+STMShield_Ttop_Pb2+0.002)/2.0,0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_left_4;
		stmShielding_left_4.name="stmShielding_left_4";
		stmShielding_left_4.solid=STMShield_left_Poly;
		finishNesting(stmShielding_left_4,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_left+G4ThreeVector((STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+STMShield_Ttop_Poly+0.003)/2.0,0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_left_5;
		stmShielding_left_5.name="stmShielding_left_5";
		stmShielding_left_5.solid=STMShield_left_Pb2;
		finishNesting(stmShielding_left_5,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_left+G4ThreeVector((STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+0.004)/2.0,0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_left_6;
		stmShielding_left_6.name="stmShielding_left_6";
		stmShielding_left_6.solid=STMShield_left_Poly;
		finishNesting(stmShielding_left_6,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_left+G4ThreeVector((STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly+2*STMShield_Ttop_Pb2+STMShield_Ttop_Poly+0.005)/2.0,0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
//for right
G4Para *STMShield_right_Al=new G4Para("STMShield_right_Al",0.5*STMShield_Tbot_Al/cos(STMShield_right_Theta),0.5*STMShield_Y1-0.001,0.5*STMShield_Z2,0,-STMShield_right_Theta,0);
G4Para *STMShield_right_Cu=new G4Para("STMShield_right_Cu",0.5*STMShield_Ttop_Cu/cos(STMShield_right_Theta),0.5*STMShield_Y1-0.001,0.5*STMShield_Z2,0,-STMShield_right_Theta,0);
G4Para *STMShield_right_Pb2=new G4Para("STMShield_right_Pb2",0.5*STMShield_Ttop_Pb2/cos(STMShield_right_Theta),0.5*STMShield_Y1-0.001,0.5*STMShield_Z2,0,-STMShield_right_Theta,0);
G4Para *STMShield_right_Poly=new G4Para("STMShield_right_Poly",0.5*STMShield_Ttop_Poly/cos(STMShield_right_Theta),0.5*STMShield_Y1-0.001,0.5*STMShield_Z2,0,-STMShield_right_Theta,0);

double STMShield_Tright_CenterX=(STMShield_X1+STMShield_X2)/4.0+STMShield_Ttop_CenterX-STMShield_thick_bot/cos(STMShield_right_Theta)+STMShield_Tbot_Al/2.0/cos(STMShield_right_Theta);
double STMShield_right_z=STMShield_Z1/2.0-(STMShield_Z1-STMShield_Zgap-STMShield_thick_front)/2.0-STMShield_Zgap;
G4ThreeVector stmShieldingPositionInParentShiftX_right=stmShieldingPositionInParent+G4ThreeVector(-STMShield_Tright_CenterX+stmSSCollxshift,stmSSCollyshift,0.0)+G4ThreeVector(-STMShield_right_z*tan(STMShield_right_Theta),0.0,STMShield_right_z);//wyq
    VolumeInfo stmShielding_right_1;
		stmShielding_right_1.name="stmShielding_right_1";
		stmShielding_right_1.solid=STMShield_right_Al;
		finishNesting(stmShielding_right_1,findMaterialOrThrow("G4_Al"),0,stmShieldingPositionInParentShiftX_right,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_right_2;
		stmShielding_right_2.name="stmShielding_right_2";
		stmShielding_right_2.solid=STMShield_right_Cu;
		finishNesting(stmShielding_right_2,findMaterialOrThrow("G4_Cu"),0,stmShieldingPositionInParentShiftX_right+G4ThreeVector(-(STMShield_Tbot_Al+STMShield_Ttop_Cu+0.001)/2.0/cos(STMShield_right_Theta),0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_right_3;
		stmShielding_right_3.name="stmShielding_right_3";
		stmShielding_right_3.solid=STMShield_right_Pb2;
		finishNesting(stmShielding_right_3,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_right+G4ThreeVector(-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+STMShield_Ttop_Pb2+0.002)/2.0/cos(STMShield_right_Theta),0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_right_4;
		stmShielding_right_4.name="stmShielding_right_4";
		stmShielding_right_4.solid=STMShield_right_Poly;
		finishNesting(stmShielding_right_4,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_right+G4ThreeVector(-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+STMShield_Ttop_Poly+0.003)/2.0/cos(STMShield_right_Theta),0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_right_5;
		stmShielding_right_5.name="stmShielding_right_5";
		stmShielding_right_5.solid=STMShield_right_Pb2;
		finishNesting(stmShielding_right_5,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_right+G4ThreeVector(-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+0.004)/2.0/cos(STMShield_right_Theta),0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
    VolumeInfo stmShielding_right_6;
		stmShielding_right_6.name="stmShielding_right_6";
		stmShielding_right_6.solid=STMShield_right_Poly;
		finishNesting(stmShielding_right_6,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_right+G4ThreeVector(-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly+2*STMShield_Ttop_Pb2+STMShield_Ttop_Poly+0.005)/2.0/cos(STMShield_right_Theta),0,0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
//<
		//>for front
///		G4Box* STMShield_front_Al=new G4Box("STMShield_front_Al",0.5*(STMShield_X1+91.44),0.5*STMShield_Y1,0.5*STMShield_Tbot_Al);
///		G4Box* STMShield_front_Cu=new G4Box("STMShield_front_Cu",0.5*(STMShield_X1+91.44),0.5*STMShield_Y1,0.5*STMShield_Ttop_Cu);
///		G4Box* STMShield_front_Pb=new G4Box("STMShield_front_Pb",0.5*(STMShield_X1+91.44),0.5*STMShield_Y1,0.5*STMShield_Ttop_Pb2);
///		G4Box* STMShield_front_Poly=new G4Box("STMShield_front_Poly",0.5*(STMShield_X1+91.44),0.5*STMShield_Y1,0.5*STMShield_Ttop_Poly);
///		double STMShield_Tfront_CenterX=-0.5*((2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Cu+STMShield_Tbot_Al)*tan(fabs(STMShield_right_Theta)))+stmSSCollxshift;
///		double STMShield_Tfront_CenterZ=-STMShield_Z1/2.0+(2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Cu+0.5*STMShield_Tbot_Al);
///		G4ThreeVector stmShieldingPositionInParentShiftX_front=stmShieldingPositionInParent+G4ThreeVector(STMShield_Tfront_CenterX,0.0,STMShield_Tfront_CenterZ);//wyq
///    VolumeInfo stmShielding_front_Poly_First;
///		stmShielding_front_Poly_First.name="stmShielding_front_1";
///		stmShielding_front_Poly_First.solid=STMShield_front_Al;
///		finishNesting(stmShielding_front_Poly_First,findMaterialOrThrow("G4_Al"),0,stmShieldingPositionInParentShiftX_front,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
///    VolumeInfo stmShielding_front_Pb_Second;
///		stmShielding_front_Pb_Second.name="stmShielding_front_Pb_Second";
///		stmShielding_front_Pb_Second.solid=STMShield_front_Cu;
///		finishNesting(stmShielding_front_Pb_Second,findMaterialOrThrow("G4_Cu"),0,stmShieldingPositionInParentShiftX_front+G4ThreeVector(0,0,-(STMShield_Tbot_Al+STMShield_Ttop_Cu+0.001)/2.0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
///    VolumeInfo stmShielding_front_Poly_Third;
///		stmShielding_front_Poly_Third.name="stmShielding_front_Poly_Third";
///		stmShielding_front_Poly_Third.solid=STMShield_front_Poly;
///		finishNesting(stmShielding_front_Poly_Third,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_front+G4ThreeVector(0,0,-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+STMShield_Ttop_Poly+0.002)/2.0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
///    VolumeInfo stmShielding_front_Cu_Fourth;
///		stmShielding_front_Cu_Fourth.name="stmShielding_front_Cu_Fourth";
///		stmShielding_front_Cu_Fourth.solid=STMShield_front_Pb;
///		finishNesting(stmShielding_front_Cu_Fourth,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShiftX_front+G4ThreeVector(0,0,-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Poly+STMShield_Ttop_Pb2+0.003)/2.0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
///    VolumeInfo stmShielding_front_Al_Fifth;
///		stmShielding_front_Al_Fifth.name="stmShielding_front_Al_Fifth";
///		stmShielding_front_Al_Fifth.solid=STMShield_front_Poly;
///		finishNesting(stmShielding_front_Al_Fifth,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShiftX_front+G4ThreeVector(0,0,-(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Poly+2*STMShield_Ttop_Pb2+STMShield_Ttop_Poly+0.004)/2.0),parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
		//two tubes used to make a hole on the front shielding
		G4Tubs* hole_LaBr3_front=new G4Tubs("hole_LaBr3_front",0.0,1.5*pSTMSSCollimatorParams.hole1RadiusDnStr(),stmSSCollHalfLength1,0.0,CLHEP::twopi);
		G4Tubs* hole_HPGe_front=new G4Tubs("hole_HPGe_front",0.0,1.5*pSTMSSCollimatorParams.hole2RadiusDnStr(),stmSSCollHalfLength1,0.0,CLHEP::twopi);
		//
////		double STMShield_front_X2=STMShield_X1+STMShield_Ttop_Poly*tan(STMShield_right_Theta);
////		double STMShield_Theta_front1=-atan((STMShield_front_X2-STMShield_X1)/2.0/STMShield_Ttop_Poly);
////		G4Trap* temp_STMShield_front_Poly_First=new G4Trap("temp_STMShield_front_Poly_First",0.5*STMShield_Ttop_Poly,STMShield_Theta_front1,0,0.5*STMShield_Y1,0.5*STMShield_X1,0.5*STMShield_X1,0,0.5*STMShield_Y1,0.5*STMShield_front_X2,0.5*STMShield_front_X2,0);
////		G4SubtractionSolid* STMShield_front_Poly_First_=new G4SubtractionSolid("STMShield_front_Poly_First_",temp_STMShield_front_Poly_First,hole_LaBr3_front,0,G4ThreeVector((STMShield_X1+STMShield_front_X2)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));
////		G4SubtractionSolid* STMShield_front_Poly_First=new G4SubtractionSolid("STMShield_front_Poly_First_",STMShield_front_Poly_First_,hole_HPGe_front,0,G4ThreeVector((STMShield_X1+STMShield_front_X2)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1)-xoffset_hole1+xoffset_hole2,-stmSSCollyshift,0.0));
////
////		double STMShield_front_X3=STMShield_X1+(STMShield_Ttop_Poly+STMShield_Ttop_Pb2)*tan(STMShield_right_Theta);
////		double STMShield_Theta_front2=-atan((STMShield_front_X3-STMShield_front_X2)/2.0/STMShield_Ttop_Pb2);
////		G4Trap* temp_STMShield_front_Pb=new G4Trap("temp_STMShield_front_Pb",0.5*STMShield_Ttop_Pb2,STMShield_Theta_front2,0,0.5*STMShield_Y1,0.5*STMShield_front_X2,0.5*STMShield_front_X2,0,0.5*STMShield_Y1,0.5*STMShield_front_X3,0.5*STMShield_front_X3,0);
////		G4SubtractionSolid* STMShield_front_Pb_=new G4SubtractionSolid("STMShield_front_Pb_",temp_STMShield_front_Pb,hole_LaBr3_front,0,G4ThreeVector((STMShield_front_X2+STMShield_front_X3)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));
////		G4SubtractionSolid* STMShield_front_Pb=new G4SubtractionSolid("STMShield_front_Pb",STMShield_front_Pb_,hole_HPGe_front,0,G4ThreeVector((STMShield_front_X2+STMShield_front_X3)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1)-xoffset_hole1+xoffset_hole2,-stmSSCollyshift,0.0));
////
////		double STMShield_front_X4=STMShield_X1+(STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Poly)*tan(STMShield_right_Theta);
////		double STMShield_Theta_front3=-atan((STMShield_front_X4-STMShield_front_X3)/2.0/STMShield_Ttop_Poly);
////		G4Trap* temp_STMShield_front_Poly_Third=new G4Trap("temp_STMShield_front_Poly_Third",0.5*STMShield_Ttop_Poly,STMShield_Theta_front3,0,0.5*STMShield_Y1,0.5*STMShield_front_X3,0.5*STMShield_front_X3,0,0.5*STMShield_Y1,0.5*STMShield_front_X4,0.5*STMShield_front_X4,0);
////		G4SubtractionSolid* STMShield_front_Poly_Third_=new G4SubtractionSolid("STMShield_front_Poly_Third_",temp_STMShield_front_Poly_Third,hole_LaBr3_front,0,G4ThreeVector((STMShield_front_X4+STMShield_front_X3)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));
////		G4SubtractionSolid* STMShield_front_Poly_Third=new G4SubtractionSolid("STMShield_front_Poly_Third",STMShield_front_Poly_Third_,hole_HPGe_front,0,G4ThreeVector((STMShield_front_X4+STMShield_front_X3)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1)-xoffset_hole1+xoffset_hole2,-stmSSCollyshift,0.0));
////
////		double STMShield_front_X5=STMShield_X1+(STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Poly+STMShield_Ttop_Cu)*tan(STMShield_right_Theta);
////		double STMShield_Theta_front4=-atan((STMShield_front_X5-STMShield_front_X4)/2.0/STMShield_Ttop_Cu);
////		G4Trap* temp_STMShield_front_Cu=new G4Trap("temp_STMShield_front_Cu",0.5*STMShield_Ttop_Cu,STMShield_Theta_front4,0,0.5*STMShield_Y1,0.5*STMShield_front_X4,0.5*STMShield_front_X4,0,0.5*STMShield_Y1,0.5*STMShield_front_X5,0.5*STMShield_front_X5,0);
////		G4SubtractionSolid* STMShield_front_Cu_=new G4SubtractionSolid("STMShield_front_Cu_",temp_STMShield_front_Cu,hole_LaBr3_front,0,G4ThreeVector((STMShield_front_X5+STMShield_front_X4)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));
////		G4SubtractionSolid* STMShield_front_Cu=new G4SubtractionSolid("STMShield_front_Cu",STMShield_front_Cu_,hole_HPGe_front,0,G4ThreeVector((STMShield_front_X5+STMShield_front_X4)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1)-xoffset_hole1+xoffset_hole2,-stmSSCollyshift,0.0));
////		double STMShield_front_X6=STMShield_X1+(STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Poly+STMShield_Ttop_Cu+STMShield_Tbot_Al)*tan(STMShield_right_Theta);
////		double STMShield_Theta_front5=-atan((STMShield_front_X6-STMShield_front_X5)/2.0/STMShield_Tbot_Al);
////		G4Trap* temp_STMShield_front_Al=new G4Trap("temp_STMShield_front_Al",0.5*STMShield_Tbot_Al,STMShield_Theta_front5,0,0.5*STMShield_Y1,0.5*STMShield_front_X5,0.5*STMShield_front_X5,0,0.5*STMShield_Y1,0.5*STMShield_front_X6,0.5*STMShield_front_X6,0);
////		G4SubtractionSolid* STMShield_front_Al_=new G4SubtractionSolid("STMShield_front_Al_",temp_STMShield_front_Al,hole_LaBr3_front,0,G4ThreeVector((STMShield_front_X5+STMShield_front_X6)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));
////		G4SubtractionSolid* STMShield_front_Al=new G4SubtractionSolid("STMShield_front_Al",STMShield_front_Al_,hole_HPGe_front,0,G4ThreeVector((STMShield_front_X5+STMShield_front_X6)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1)-xoffset_hole1+xoffset_hole2,-stmSSCollyshift,0.0));
		G4UnionSolid* holes_front=new G4UnionSolid("holes_front",hole_LaBr3_front,hole_HPGe_front,0,G4ThreeVector(-xoffset_hole1+xoffset_hole2,0.0,0.0));
		//
		double STMShield_front_X2=STMShield_X1+STMShield_Ttop_Poly*tan(STMShield_right_Theta);
		double STMShield_Theta_front1=-atan((STMShield_front_X2-STMShield_X1)/2.0/STMShield_Ttop_Poly);
		G4Trap* temp_STMShield_front_Poly_First=new G4Trap("temp_STMShield_front_Poly_First",0.5*STMShield_Ttop_Poly,STMShield_Theta_front1,0,0.5*STMShield_Y1,0.5*STMShield_X1,0.5*STMShield_X1,0,0.5*STMShield_Y1,0.5*STMShield_front_X2,0.5*STMShield_front_X2,0);
		G4SubtractionSolid* STMShield_front_Poly_First=new G4SubtractionSolid("STMShield_front_Poly_First",temp_STMShield_front_Poly_First,holes_front,0,G4ThreeVector((STMShield_X1+STMShield_front_X2)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));

		double STMShield_front_X3=STMShield_X1+(STMShield_Ttop_Poly+STMShield_Ttop_Pb2)*tan(STMShield_right_Theta);
		double STMShield_Theta_front2=-atan((STMShield_front_X3-STMShield_front_X2)/2.0/STMShield_Ttop_Pb2);
		G4Trap* temp_STMShield_front_Pb=new G4Trap("temp_STMShield_front_Pb",0.5*STMShield_Ttop_Pb2,STMShield_Theta_front2,0,0.5*STMShield_Y1,0.5*STMShield_front_X2,0.5*STMShield_front_X2,0,0.5*STMShield_Y1,0.5*STMShield_front_X3,0.5*STMShield_front_X3,0);
		G4SubtractionSolid* STMShield_front_Pb=new G4SubtractionSolid("STMShield_front_Pb",temp_STMShield_front_Pb,holes_front,0,G4ThreeVector((STMShield_front_X2+STMShield_front_X3)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));

		double STMShield_front_X4=STMShield_X1+(STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Poly)*tan(STMShield_right_Theta);
		double STMShield_Theta_front3=-atan((STMShield_front_X4-STMShield_front_X3)/2.0/STMShield_Ttop_Poly);
		G4Trap* temp_STMShield_front_Poly_Third=new G4Trap("temp_STMShield_front_Poly_Third",0.5*STMShield_Ttop_Poly,STMShield_Theta_front3,0,0.5*STMShield_Y1,0.5*STMShield_front_X3,0.5*STMShield_front_X3,0,0.5*STMShield_Y1,0.5*STMShield_front_X4,0.5*STMShield_front_X4,0);
		G4SubtractionSolid* STMShield_front_Poly_Third=new G4SubtractionSolid("STMShield_front_Poly_Third",temp_STMShield_front_Poly_Third,holes_front,0,G4ThreeVector((STMShield_front_X4+STMShield_front_X3)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));

		double STMShield_front_X5=STMShield_X1+(STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Poly+STMShield_Ttop_Cu)*tan(STMShield_right_Theta);
		double STMShield_Theta_front4=-atan((STMShield_front_X5-STMShield_front_X4)/2.0/STMShield_Ttop_Cu);
		G4Trap* temp_STMShield_front_Cu=new G4Trap("temp_STMShield_front_Cu",0.5*STMShield_Ttop_Cu,STMShield_Theta_front4,0,0.5*STMShield_Y1,0.5*STMShield_front_X4,0.5*STMShield_front_X4,0,0.5*STMShield_Y1,0.5*STMShield_front_X5,0.5*STMShield_front_X5,0);
		G4SubtractionSolid* STMShield_front_Cu=new G4SubtractionSolid("STMShield_front_Cu",temp_STMShield_front_Cu,holes_front,0,G4ThreeVector((STMShield_front_X5+STMShield_front_X4)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));
		double STMShield_front_X6=STMShield_X1+(STMShield_Ttop_Poly+STMShield_Ttop_Pb2+STMShield_Ttop_Poly+STMShield_Ttop_Cu+STMShield_Tbot_Al)*tan(STMShield_right_Theta);
		double STMShield_Theta_front5=-atan((STMShield_front_X6-STMShield_front_X5)/2.0/STMShield_Tbot_Al);
		G4Trap* temp_STMShield_front_Al=new G4Trap("temp_STMShield_front_Al",0.5*STMShield_Tbot_Al,STMShield_Theta_front5,0,0.5*STMShield_Y1,0.5*STMShield_front_X5,0.5*STMShield_front_X5,0,0.5*STMShield_Y1,0.5*STMShield_front_X6,0.5*STMShield_front_X6,0);
		G4SubtractionSolid* STMShield_front_Al=new G4SubtractionSolid("STMShield_front_Al",temp_STMShield_front_Al,holes_front,0,G4ThreeVector((STMShield_front_X5+STMShield_front_X6)/4.0-(stmSSCollHalfWidth1+stmSSCollxshift-xoffset_hole1),-stmSSCollyshift,0.0));
		double STMShield_front_xshift=STMShield_Ttop_Poly*0.5*tan(STMShield_Theta_front1);
		double STMShield_front_zshift=-STMShield_Z1/2.0+0.5*STMShield_Ttop_Poly;
		G4ThreeVector stmShieldingPositionInParentShift_front=stmShieldingPositionInParent+G4ThreeVector(STMShield_front_xshift,0.0,STMShield_front_zshift)+stmSSCollShift;
    VolumeInfo stmShielding_front_Poly_First;
		stmShielding_front_Poly_First.name="stmShielding_front_Poly_First";
		stmShielding_front_Poly_First.solid=STMShield_front_Poly_First;
		finishNesting(stmShielding_front_Poly_First,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShift_front,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);

		STMShield_front_xshift=STMShield_Ttop_Pb2*0.5*tan(STMShield_Theta_front2)-(STMShield_front_X2-STMShield_X1)/2.0;
		STMShield_front_zshift=-STMShield_Z1/2.0+0.5*(STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly);
		stmShieldingPositionInParentShift_front=stmShieldingPositionInParent+G4ThreeVector(STMShield_front_xshift,0.0,STMShield_front_zshift)+stmSSCollShift;
    VolumeInfo stmShielding_front_Pb_Second;
		stmShielding_front_Pb_Second.name="stmShielding_front_Pb_Second";
		stmShielding_front_Pb_Second.solid=STMShield_front_Pb;
		finishNesting(stmShielding_front_Pb_Second,findMaterialOrThrow("G4_Pb"),0,stmShieldingPositionInParentShift_front,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);

		STMShield_front_xshift=STMShield_Ttop_Poly*0.5*tan(STMShield_Theta_front3)-(STMShield_front_X3-STMShield_X1)/2.0;
		STMShield_front_zshift=-STMShield_Z1/2.0+0.5*(STMShield_Ttop_Poly+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly);
		stmShieldingPositionInParentShift_front=stmShieldingPositionInParent+G4ThreeVector(STMShield_front_xshift,0.0,STMShield_front_zshift)+stmSSCollShift;
    VolumeInfo stmShielding_front_Poly_Third;
		stmShielding_front_Poly_Third.name="stmShielding_front_Poly_Third";
		stmShielding_front_Poly_Third.solid=STMShield_front_Poly_Third;
		finishNesting(stmShielding_front_Poly_Third,findMaterialOrThrow("Polyethylene092B050d095"),0,stmShieldingPositionInParentShift_front,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);

		STMShield_front_xshift=STMShield_Ttop_Cu*0.5*tan(STMShield_Theta_front4)-(STMShield_front_X4-STMShield_X1)/2.0;
		STMShield_front_zshift=-STMShield_Z1/2.0+0.5*(STMShield_Ttop_Cu+2*STMShield_Ttop_Poly+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly);
		stmShieldingPositionInParentShift_front=stmShieldingPositionInParent+G4ThreeVector(STMShield_front_xshift,0.0,STMShield_front_zshift)+stmSSCollShift;

    VolumeInfo stmShielding_front_Cu_Fourth;
		stmShielding_front_Cu_Fourth.name="stmShielding_front_Cu_Fourth";
		stmShielding_front_Cu_Fourth.solid=STMShield_front_Cu;
		finishNesting(stmShielding_front_Cu_Fourth,findMaterialOrThrow("G4_Cu"),0,stmShieldingPositionInParentShift_front,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);

		STMShield_front_xshift=STMShield_Tbot_Al*0.5*tan(STMShield_Theta_front5)-(STMShield_front_X5-STMShield_X1)/2.0;
		STMShield_front_zshift=-STMShield_Z1/2.0+0.5*(STMShield_Tbot_Al+2*STMShield_Ttop_Cu+2*STMShield_Ttop_Poly+2*STMShield_Ttop_Pb2+2*STMShield_Ttop_Poly);
		stmShieldingPositionInParentShift_front=stmShieldingPositionInParent+G4ThreeVector(STMShield_front_xshift,0.0,STMShield_front_zshift)+stmSSCollShift;
    VolumeInfo stmShielding_front_Al_Fifth;
		stmShielding_front_Al_Fifth.name="stmShielding_front_Al_Fifth";
		stmShielding_front_Al_Fifth.solid=STMShield_front_Al;
		finishNesting(stmShielding_front_Al_Fifth,findMaterialOrThrow("G4_Al"),0,stmShieldingPositionInParentShift_front,parentInfo.logical,0,STMisVisible,G4Color::Cyan(),STMisSolid,forceAuxEdgeVisible,placePV,doSurfaceCheck);
		//<
    G4Material*  stmDet2Material                 =  findMaterialOrThrow(pSTMDetector2Params.crystalMaterial());
    const double stmDet2ROut                     =  pSTMDetector2Params.crystalRadiusOut();
    const double stmDet2HalfLength               =  pSTMDetector2Params.crystalHalfLength();
///    G4Material*  stmDet2CanMaterial              =  findMaterialOrThrow(pSTMDetector2Params.canMaterial());
    const double stmDet2CanRIn                   =  pSTMDetector2Params.canRadiusIn();
    const double stmDet2CanROut                  =  pSTMDetector2Params.canRadiusOut();
    const double stmDet2CanHalfLength            =  pSTMDetector2Params.canHalfLength();
    const double stmDet2CanUpStrWindowHalfLength =  pSTMDetector2Params.canUpStrWindowHalfLength();
///    G4Material*  stmDet2CanDnStrWindowMaterial   =  stmDet2CanMaterial;
    const double stmDet2CanDnStrWindowHalfLength =  stmDet2CanROut - stmDet2CanRIn;
    const TubsParams stmDet2Params(0., stmDet2ROut, stmDet2HalfLength);
    const TubsParams stmDet2CanParams(stmDet2CanRIn, stmDet2CanROut, stmDet2CanHalfLength);
    const TubsParams stmDet2CanGasParams(0.,  stmDet2CanRIn, stmDet2CanHalfLength - stmDet2CanUpStrWindowHalfLength - stmDet2CanDnStrWindowHalfLength);
    const TubsParams stmDet2CanUpStrWindowParams(0., stmDet2CanRIn, stmDet2CanUpStrWindowHalfLength);
    const TubsParams stmDet2CanDnStrWindowParams(0., stmDet2CanRIn, stmDet2CanDnStrWindowHalfLength);

    G4ThreeVector stmHPGePositionInMu2e   = pSTMDetector2Params.originInMu2e();
    G4ThreeVector stmHPGePositionInParent = pSTMDetector2Params.originInMu2e() - parentCenterInMu2e;

    G4Tubs *stmDet2outPart1=new G4Tubs("stmDet2outPart1",0.0,STMDet2HPGe_R_out,0.5*(STMDet2HPGe_L_out-STMDet2HPGe_R_out_corner),0.0,CLHEP::twopi );//72.1/2 for radius, 78.5-8 for length
		G4Torus *stmDet2outPart2=new G4Torus("stmDet2outPart2",0.0,STMDet2HPGe_R_out_corner,STMDet2HPGe_R_out-STMDet2HPGe_R_out_corner,0.0,CLHEP::twopi);//36.05-8
    G4Tubs *stmDet2outPart3=new G4Tubs("stmDet2outPart3",0.0,STMDet2HPGe_R_out-STMDet2HPGe_R_out_corner,STMDet2HPGe_R_out_corner,0.0,CLHEP::twopi);
    G4UnionSolid *stmDet2outPart2_3=new G4UnionSolid("stmDet2outPart2_3",stmDet2outPart2,stmDet2outPart3,0,G4ThreeVector(0.0,0.0,0.0));
    G4UnionSolid *stmDet2outPart1_2_3=new G4UnionSolid("stmDet2outPart1_2_3",stmDet2outPart1,stmDet2outPart2_3,0,G4ThreeVector(0.0,0.0,-0.5*(STMDet2HPGe_L_out-STMDet2HPGe_R_out_corner)));
    G4Tubs *stmDet2inPart1=new G4Tubs("stmDet2inPart1",0.0,STMDet2HPGe_R_in,0.5*(STMDet2HPGe_L_in-STMDet2HPGe_R_in_corner),0.0,CLHEP::twopi );//72.1/2 for radius, 78.5-8 for length
		G4Sphere *stmDet2inPart2=new G4Sphere("stmDet2inPart2",0.0,STMDet2HPGe_R_in_corner,0.0,CLHEP::twopi,0.0,CLHEP::twopi);//36.05-8
    G4UnionSolid *stmDet2inPart1_2=new G4UnionSolid("stmDet2inPart1_2",stmDet2inPart1,stmDet2inPart2,0,G4ThreeVector(0.0,0.0,-0.5*(STMDet2HPGe_L_in-STMDet2HPGe_R_in_corner)));
    G4SubtractionSolid *STMDet2HPGe=new G4SubtractionSolid("STMDet2HPGe",stmDet2outPart1_2_3,stmDet2inPart1_2,0,G4ThreeVector(0.0,0.0,0.5*(STMDet2HPGe_L_out-STMDet2HPGe_R_out_corner-STMDet2HPGe_L_in+STMDet2HPGe_R_in_corner)+0.001));
//>cup for HPGe
//double STMDet2HPGe_Cup0_gap_top=0.0;
double STMDet2HPGe_Cup0_halfLength=(STMDet2HPGe_L_out+STMDet2HPGe_Cup0_thick_top+STMDet2HPGe_Cup0_thick_bot+STMDet2HPGe_Cup0_gap_bot)/2.0;
double STMDet2HPGe_Cup0_R=STMDet2HPGe_R_out+STMDet2HPGe_Cup0_gap_bar+STMDet2HPGe_Cup0_thick_bar;
G4Tubs* STMDet2HPGe_Cup0=new G4Tubs("STMDet2HPGe_Cup0",0,STMDet2HPGe_Cup0_R,STMDet2HPGe_Cup0_halfLength,0.0,CLHEP::twopi);
G4Tubs* STMDet2HPGe_Cup1=new G4Tubs("STMDet2HPGe_Cup1",0,STMDet2HPGe_Cup0_R-STMDet2HPGe_Cup0_thick_bar,STMDet2HPGe_Cup0_halfLength-(STMDet2HPGe_Cup0_thick_top+STMDet2HPGe_Cup0_thick_bot)/2.0,0.0,CLHEP::twopi);
G4SubtractionSolid* STMDet2HPGe_Cup2=new G4SubtractionSolid("STMDet2HPGe_Cup2",STMDet2HPGe_Cup0,STMDet2HPGe_Cup1,0,G4ThreeVector(0.0,0.0,(-STMDet2HPGe_Cup0_thick_bot+STMDet2HPGe_Cup0_thick_top)/2.0));
G4SubtractionSolid* STMDet2HPGe_Cup=new G4SubtractionSolid("STMDet2HPGe_Cup",STMDet2HPGe_Cup2,stmDet2inPart1,0,G4ThreeVector(0.0,0.0,STMDet2HPGe_L_in));
//>Al can
double STMDet2HPGe_Can_gap_bot=STMDet2HPGe_Can_halfL*2-2*STMDet2HPGe_Cup0_halfLength-STMDet2HPGe_Can_gap_top-STMDet2HPGe_Can_thick;
G4Tubs* STMDet2HPGe_Can0=new G4Tubs("STMDet2HPGe_Can0",0,STMDet2HPGe_Can_R,STMDet2HPGe_Can_halfL,0.0,CLHEP::twopi);
G4Tubs* STMDet2HPGe_Can1=new G4Tubs("STMDet2HPGe_Can1",0,STMDet2HPGe_Can_R-STMDet2HPGe_Can_thick,STMDet2HPGe_Can_halfL,0.0,CLHEP::twopi);
G4SubtractionSolid* STMDet2HPGe_Can=new G4SubtractionSolid("STMDet2HPGe_Can",STMDet2HPGe_Can0,STMDet2HPGe_Can1,0,G4ThreeVector(0.0,0.0,STMDet2HPGe_Can_thick));
//>Positions for the HPGe related
CLHEP::HepRotationY angleHPGe(STMShield_right_Theta);
G4RotationMatrix *rotHPGe = new G4RotationMatrix(angleHPGe);
G4ThreeVector stmDet2CanPositionInParent=stmHPGePositionInParent+G4ThreeVector(STMDet2HPGe_Can_halfL*sin(-STMShield_right_Theta)+STMDet2HPGe_xOffset,0.0,-STMDet2HPGe_Can_halfL*(1-cos(STMShield_right_Theta))+STMShield_thick_front);//+stmSSCollShift;
    VolumeInfo stmDet2Can;
    stmDet2Can.name="stmDet2Can";
    stmDet2Can.solid=STMDet2HPGe_Can;
    finishNesting(stmDet2Can,
				findMaterialOrThrow("G4_Al"),rotHPGe,stmDet2CanPositionInParent,
       parentInfo.logical, 0, STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);

double shiftDistanceXY_Cup=STMDet2HPGe_Can_gap_bot/2.0-(STMDet2HPGe_Can_thick+STMDet2HPGe_Can_gap_top)/2.0;
G4ThreeVector stmDet2CupPositionInParent=stmDet2CanPositionInParent+G4ThreeVector(shiftDistanceXY_Cup*sin(STMShield_right_Theta),0.0,-shiftDistanceXY_Cup*cos(STMShield_right_Theta));
    VolumeInfo stmDet2Cup;
    stmDet2Cup.name="stmDet2Cup";
    stmDet2Cup.solid=STMDet2HPGe_Cup;
    finishNesting(stmDet2Cup,
				findMaterialOrThrow("G4_Al"),rotHPGe,stmDet2CupPositionInParent,
       parentInfo.logical, 0, STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);

double shiftDistanceXY_Ge=(STMDet2HPGe_Can_gap_bot+STMDet2HPGe_Cup0_thick_bot+STMDet2HPGe_Cup0_gap_bot)/2.0-(STMDet2HPGe_Can_thick+STMDet2HPGe_Can_gap_top+STMDet2HPGe_Cup0_thick_top+STMDet2HPGe_R_out_corner+0.001)/2.0;
G4ThreeVector stmDet2PositionInParent=stmDet2CanPositionInParent+G4ThreeVector(shiftDistanceXY_Ge*sin(STMShield_right_Theta),0.0,-shiftDistanceXY_Ge*cos(STMShield_right_Theta));
    VolumeInfo stmDet2;
    stmDet2.name="stmDet2";
    stmDet2.solid=STMDet2HPGe;
    finishNesting(stmDet2,stmDet2Material,rotHPGe,
			 stmDet2PositionInParent, parentInfo.logical, 1,
       STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
    G4VSensitiveDetector* sd2=G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::STMDet());
  if(sd2)stmDet2.logical->SetSensitiveDetector(sd2);
//define inner solids
G4Tubs* STMShield_LaBr3_hole=new G4Tubs("STMShield_LaBr3_hole",0.0,STMShield_LaBr3_Hole_R,STMShield_LaBr3_halfT_Pb*1.1,0.0,CLHEP::twopi);//make the hole
G4Box* _STMShield_LaBr3_Al=new G4Box("_STMShield_LaBr3_Al",STMShield_LaBr3_halfWidth,STMShield_Y1/2.0,STMShield_LaBr3_halfT_Al);
G4Box* _STMShield_LaBr3_Cu=new G4Box("_STMShield_LaBr3_Cu",STMShield_LaBr3_halfWidth,STMShield_Y1/2.0,STMShield_LaBr3_halfT_Cu);
G4Box* _STMShield_LaBr3_Pb=new G4Box("_STMShield_LaBr3_Pb",STMShield_LaBr3_halfWidth,STMShield_Y1/2.0,STMShield_LaBr3_halfT_Pb);
G4SubtractionSolid* STMShield_LaBr3_Al=new G4SubtractionSolid("STMShield_LaBr3_Al",_STMShield_LaBr3_Al,STMShield_LaBr3_hole,0,G4ThreeVector(STMShield_thick_bot+STMShield_LaBr3_halfWidth-(stmSSCollHalfWidth1-xoffset_hole1+stmSSCollxshift),-stmSSCollyshift,0.0));
G4SubtractionSolid* STMShield_LaBr3_Cu=new G4SubtractionSolid("STMShield_LaBr3_Cu",_STMShield_LaBr3_Cu,STMShield_LaBr3_hole,0,G4ThreeVector(STMShield_thick_bot+STMShield_LaBr3_halfWidth-(stmSSCollHalfWidth1-xoffset_hole1+stmSSCollxshift),-stmSSCollyshift,0.0));
G4SubtractionSolid* STMShield_LaBr3_Pb=new G4SubtractionSolid("STMShield_LaBr3_Pb",_STMShield_LaBr3_Pb,STMShield_LaBr3_hole,0,G4ThreeVector(STMShield_thick_bot+STMShield_LaBr3_halfWidth-(stmSSCollHalfWidth1-xoffset_hole1+stmSSCollxshift),-stmSSCollyshift,0.0));

//start locations
G4ThreeVector stmShield_LaBr3_AlPositionInParent=stmShieldingPositionInParentShiftX_left+G4ThreeVector(-STMShield_Tbot_Al/2.0-STMShield_LaBr3_halfWidth,0.0,-0.5*STMShield_Z2+STMShield_LaBr3_gap_Al+STMShield_LaBr3_halfT_Al);
    VolumeInfo stmShield_LaBr3_Al;
    stmShield_LaBr3_Al.name="stmShield_LaBr3_Al";
    stmShield_LaBr3_Al.solid=STMShield_LaBr3_Al;
    finishNesting(stmShield_LaBr3_Al,findMaterialOrThrow("G4_Al"),0,
				stmShield_LaBr3_AlPositionInParent,parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);

    VolumeInfo stmShield_LaBr3_Cu;
    stmShield_LaBr3_Cu.name="stmShield_LaBr3_Cu";
    stmShield_LaBr3_Cu.solid=STMShield_LaBr3_Cu;
    finishNesting(stmShield_LaBr3_Cu,findMaterialOrThrow("G4_Cu"),0,
				stmShield_LaBr3_AlPositionInParent+G4ThreeVector(0.0,0.0,STMShield_LaBr3_halfT_Al+STMShield_LaBr3_halfT_Cu),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);

    VolumeInfo stmShield_LaBr3_Pb;
    stmShield_LaBr3_Pb.name="stmShield_LaBr3_Pb";
    stmShield_LaBr3_Pb.solid=STMShield_LaBr3_Pb;
    finishNesting(stmShield_LaBr3_Pb,findMaterialOrThrow("G4_Pb"),0,
				stmShield_LaBr3_AlPositionInParent+G4ThreeVector(0.0,0.0,STMShield_LaBr3_halfT_Al+STMShield_LaBr3_halfT_Pb+2*STMShield_LaBr3_halfT_Cu),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);

    VolumeInfo stmShield2_LaBr3_Al;
    stmShield2_LaBr3_Al.name="stmShield2_LaBr3_Al";
    stmShield2_LaBr3_Al.solid=STMShield_LaBr3_Al;
    finishNesting(stmShield2_LaBr3_Al,findMaterialOrThrow("G4_Al"),0,
				stmShield_LaBr3_AlPositionInParent+G4ThreeVector(0.0,0.0,2*STMShield_LaBr3_halfT_Al+4*STMShield_LaBr3_halfT_Cu+4*STMShield_LaBr3_halfT_Pb+STMShield_LaBr3_gap_Pb),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);

    VolumeInfo stmShield2_LaBr3_Cu;
    stmShield2_LaBr3_Cu.name="stmShield2_LaBr3_Cu";
    stmShield2_LaBr3_Cu.solid=STMShield_LaBr3_Cu;
    finishNesting(stmShield2_LaBr3_Cu,findMaterialOrThrow("G4_Cu"),0,
				stmShield_LaBr3_AlPositionInParent+G4ThreeVector(0.0,0.0,STMShield_LaBr3_halfT_Al+STMShield_LaBr3_halfT_Cu*3+4*STMShield_LaBr3_halfT_Pb+STMShield_LaBr3_gap_Pb),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);

    VolumeInfo stmShield2_LaBr3_Pb;
    stmShield2_LaBr3_Pb.name="stmShield2_LaBr3_Pb";
    stmShield2_LaBr3_Pb.solid=STMShield_LaBr3_Pb;
    finishNesting(stmShield2_LaBr3_Pb,findMaterialOrThrow("G4_Pb"),0,
				stmShield_LaBr3_AlPositionInParent+G4ThreeVector(0.0,0.0,STMShield_LaBr3_halfT_Al+STMShield_LaBr3_halfT_Pb*3+2*STMShield_LaBr3_halfT_Cu+STMShield_LaBr3_gap_Pb),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
//>shielding btw. HPGe and LaBr3
G4Box* STMShield_I_Pb=new G4Box("STMShield_I_Pb",STMShield_halfTI_Pb,STMShield_Y1/2.0,STMShield_halfLengthI);
G4Box* STMShield_I_Al=new G4Box("STMShield_I_Al",STMShield_halfTI_Al,STMShield_Y1/2.0,STMShield_halfLengthI);
G4Box* STMShield_I_Cu=new G4Box("STMShield_I_Cu",STMShield_halfTI_Cu,STMShield_Y1/2.0,STMShield_halfLengthI);
G4ThreeVector stmShield_I_AlPositionInParent=stmShield_LaBr3_AlPositionInParent+G4ThreeVector(-STMShield_LaBr3_halfWidth-STMShield_halfTI_Al,0.0,-STMShield_LaBr3_halfT_Al+STMShield_halfLengthI);
//locations
    VolumeInfo stmShield_I_Al;
    stmShield_I_Al.name="stmShield_I_Al";
    stmShield_I_Al.solid=STMShield_I_Al;
    finishNesting(stmShield_I_Al,findMaterialOrThrow("G4_Al"),0,
				stmShield_I_AlPositionInParent,parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
    VolumeInfo stmShield_I_Cu;
    stmShield_I_Cu.name="stmShield_I_Cu";
    stmShield_I_Cu.solid=STMShield_I_Cu;
    finishNesting(stmShield_I_Cu,findMaterialOrThrow("G4_Cu"),0,
				stmShield_I_AlPositionInParent+G4ThreeVector(-STMShield_halfTI_Al-STMShield_halfTI_Cu,0.0,0.0),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
    VolumeInfo stmShield_I_Pb;
    stmShield_I_Pb.name="stmShield_I_Pb";
    stmShield_I_Pb.solid=STMShield_I_Pb;
    finishNesting(stmShield_I_Pb,findMaterialOrThrow("G4_Pb"),0,
				stmShield_I_AlPositionInParent+G4ThreeVector(-STMShield_halfTI_Al-2*STMShield_halfTI_Cu-STMShield_halfTI_Pb,0.0,0.0),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
    VolumeInfo stmShield_I_Cu2;
    stmShield_I_Cu2.name="stmShield_I_Cu2";
    stmShield_I_Cu2.solid=STMShield_I_Cu;
    finishNesting(stmShield_I_Cu2,findMaterialOrThrow("G4_Cu"),0,
				stmShield_I_AlPositionInParent+G4ThreeVector(-STMShield_halfTI_Al-STMShield_halfTI_Cu*3-STMShield_halfTI_Pb*2,0.0,0.0),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
    VolumeInfo stmShield_I_Al2;
    stmShield_I_Al2.name="stmShield_I_Al2";
    stmShield_I_Al2.solid=STMShield_I_Al;
    finishNesting(stmShield_I_Al2,findMaterialOrThrow("G4_Al"),0,
				stmShield_I_AlPositionInParent+G4ThreeVector(-2*STMShield_halfTI_Al-STMShield_halfTI_Cu*4-STMShield_halfTI_Pb*2,0.0,0.0),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
//>LaBr3
G4Tubs* _STMDet1CupLaBr3=new G4Tubs("_STMDet1CupLaBr3",0.0,STMDet1CupLaBr3_R,STMDet1CupLaBr3_halfLength,0.,CLHEP::twopi);
G4Tubs* _subSTMDet1CupLaBr3=new G4Tubs("_subSTMDet1CupLaBr3",0.0,STMDet1CupLaBr3_R-STMDet1CupLaBr3_thick,STMDet1CupLaBr3_halfLength,0.0,CLHEP::twopi);
G4SubtractionSolid* STMDet1CupLaBr3=new G4SubtractionSolid("STMDet1CupLaBr3",_STMDet1CupLaBr3,_subSTMDet1CupLaBr3,0,G4ThreeVector(0.0,0.0,STMDet1CupLaBr3_thick));
///G4ThreeVector stmDet1CupPositionInParent=G4ThreeVector(stmShield_I_AlPositionInParent.x()+STMShield_halfTI_Al+STMDet1CupLaBr3_R,stmShield_I_AlPositionInParent.y(),stmShield_LaBr3_AlPositionInParent.z()+4*(STMShield_LaBr3_halfT_Cu+STMShield_LaBr3_halfT_Pb)+STMShield_LaBr3_halfT_Al*3+STMDet1CupLaBr3_halfLength+STMShield_LaBr3_gap_Pb+STMDet1CupLaBr3_Upstream);
G4ThreeVector stmDet1CupPositionInParent=stmShieldingPositionInParent+G4ThreeVector(STMDet1LaBr3_xOffset+xoffset_hole1,0.0,STMShield_thick_front+STMShield_LaBr3_gap_Al+4*(STMShield_LaBr3_halfT_Cu+STMShield_LaBr3_halfT_Pb+STMShield_LaBr3_halfT_Al)+STMDet1CupLaBr3_halfLength+STMShield_LaBr3_gap_Pb+STMDet1CupLaBr3_Upstream-STMShield_Z1/2.0);
G4Tubs* STMDet1LaBr3=new G4Tubs("STMDet1LaBr3",0.0,STMDet1LaBr3_R,STMDet1LaBr3_halfLength,0.,CLHEP::twopi);//Crystal

    VolumeInfo stmDet1CupLaBr3_Al;
    stmDet1CupLaBr3_Al.name="stmDet1CupLaBr3_Al";
    stmDet1CupLaBr3_Al.solid=STMDet1CupLaBr3;
    finishNesting(stmDet1CupLaBr3_Al,findMaterialOrThrow("G4_Al"),0,
				stmDet1CupPositionInParent,parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
    VolumeInfo stmDet1;
    stmDet1.name="stmDet1";
    stmDet1.solid=STMDet1LaBr3;
    finishNesting(stmDet1,findMaterialOrThrow("LaBr3Ce"),0,
				stmDet1CupPositionInParent+G4ThreeVector(0.0,0.0,-STMDet1CupLaBr3_halfLength+STMDet1CupLaBr3_thick+STMDet1LaBr3_gap_Upstream+STMDet1LaBr3_halfLength),parentInfo.logical, 0,
				STMisVisible, G4Color::Gray(), STMisSolid, forceAuxEdgeVisible, placePV, doSurfaceCheck);
    // Make stmDet1 a sensitive detector.
    G4VSensitiveDetector* sd1 = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::STMDet());
    if(sd1) stmDet1.logical->SetSensitiveDetector(sd1);

///    if(pSTMDetector2Params.build()){
///    G4Material*  stmDet2Material                 =  findMaterialOrThrow(pSTMDetector2Params.crystalMaterial());
///    const double stmDet2ROut                     =  pSTMDetector2Params.crystalRadiusOut();
///    const double stmDet2HalfLength               =  pSTMDetector2Params.crystalHalfLength();
///    G4Material*  stmDet2CanMaterial              =  findMaterialOrThrow(pSTMDetector2Params.canMaterial());
///    const double stmDet2CanRIn                   =  pSTMDetector2Params.canRadiusIn();
///    const double stmDet2CanROut                  =  pSTMDetector2Params.canRadiusOut();
///    const double stmDet2CanHalfLength            =  pSTMDetector2Params.canHalfLength();
///    G4Material*  stmDet2CanUpStrWindowMaterial   =  findMaterialOrThrow(pSTMDetector2Params.canUpStrWindowMaterial());
///    const double stmDet2CanUpStrWindowHalfLength =  pSTMDetector2Params.canUpStrWindowHalfLength();
///    G4Material*  stmDet2CanDnStrWindowMaterial   =  stmDet2CanMaterial;
///    const double stmDet2CanDnStrWindowHalfLength =  stmDet2CanROut - stmDet2CanRIn;
///
///    if ( stmDet2ROut > stmDet2CanRIn ){
///      throw cet::exception("GEOM")<< " STM: det1 radius is larger than the inner radius of the can. \n" ;
///    }
///    if ( stmDet2HalfLength > stmDet2CanHalfLength-stmDet2CanDnStrWindowHalfLength-stmDet2CanUpStrWindowHalfLength ){
///      throw cet::exception("GEOM")<< " STM: det1 crystal length is larger than the inner length of the can. \n" ;
///    }
///
///    const TubsParams stmDet2Params(0., stmDet2ROut, stmDet2HalfLength);
///    const TubsParams stmDet2CanParams(stmDet2CanRIn, stmDet2CanROut, stmDet2CanHalfLength);
///    const TubsParams stmDet2CanGasParams(0.,  stmDet2CanRIn, stmDet2CanHalfLength - stmDet2CanUpStrWindowHalfLength - stmDet2CanDnStrWindowHalfLength);
///    const TubsParams stmDet2CanUpStrWindowParams(0., stmDet2CanRIn, stmDet2CanUpStrWindowHalfLength);
///    const TubsParams stmDet2CanDnStrWindowParams(0., stmDet2CanRIn, stmDet2CanDnStrWindowHalfLength);
///
///    G4ThreeVector stmDet2CanPositionInMu2e   = pSTMDetector2Params.originInMu2e();
///    G4ThreeVector stmDet2CanPositionInParent = pSTMDetector2Params.originInMu2e() - parentCenterInMu2e;
///    G4ThreeVector stmHPGePositionInMu2e      = stmDet2CanPositionInMu2e;
///    G4ThreeVector stmHPGePositionInParent    = stmDet2CanPositionInParent;
///
///    VolumeInfo stmDet2CanInfo = nestTubs( "stmDet2Can",
///                                         stmDet2CanParams,
///                                         stmDet2CanMaterial,
///                                         0x0,
///                                         stmDet2CanPositionInParent,
///                                         parentInfo,
///                                         0,
///                                         STMisVisible,
///                                         G4Color::Green(),
///                                         STMisSolid,
///                                         forceAuxEdgeVisible,
///                                         placePV,
///                                         doSurfaceCheck
///                                         );
///
///    VolumeInfo stmDet2 = nestTubs("stmDet2",
///                                      stmDet2Params,
///                                      stmDet2Material,
///                                      0x0,
///                                      stmDet2CanPositionInParent,
///                                      parentInfo,
///                                      1,
///                                      STMisVisible,
///                                      G4Color::Red(),
///                                      STMisSolid,
///                                      forceAuxEdgeVisible,
///                                      placePV,
///                                      doSurfaceCheck
///                                      );
///
///    VolumeInfo stmDet2CanUpStrWindowInfo = nestTubs( "stmDet2CanUpStrWindow",
///                                            stmDet2CanUpStrWindowParams,
///                                            stmDet2CanUpStrWindowMaterial,
///                                            0x0,
///                                            stmDet2CanPositionInParent + G4ThreeVector(0.0,0.0,-1.0*stmDet2CanHalfLength + stmDet2CanUpStrWindowHalfLength),
///                                            parentInfo,
///                                            0,
///                                            STMisVisible,
///                                            G4Color::Green(),
///                                            STMisSolid,
///                                            forceAuxEdgeVisible,
///                                            placePV,
///                                            doSurfaceCheck
///                                            );
///
///    VolumeInfo stmDet2CanDnStrWindowInfo = nestTubs( "stmDet2CanDnStrWindow",
///                                            stmDet2CanDnStrWindowParams,
///                                            stmDet2CanDnStrWindowMaterial,
///                                            0x0,
///                                            stmDet2CanPositionInParent + G4ThreeVector(0.0,0.0, stmDet2CanHalfLength - stmDet2CanDnStrWindowHalfLength),
///                                            parentInfo,
///                                            0,
///                                            STMisVisible,
///                                            G4Color::Green(),
///                                            STMisSolid,
///                                            forceAuxEdgeVisible,
///                                            placePV,
///                                            doSurfaceCheck
///                                            );
///
///    if (verbosityLevel>0){
///      std::cout << __func__ << " Warning: Gas not implemented inside STM detector1 can! (so that VD inside can does not overlap with can gas)" << std::endl;
///    }
///
///    // Make stmDet2 a sensitive detector.
///    G4VSensitiveDetector *sd2 = G4SDManager::GetSDMpointer()->FindSensitiveDetector(SensitiveDetectorName::STMDet());
///    if(sd2) stmDet2.logical->SetSensitiveDetector(sd2);
///		}

    //===================== Shield Pipe/Wall to prevent michel electrons from causing deadtime in the CRV  ==========================

    const double mstmCRVShieldDnStrSpace     =  pSTMShieldPipeParams.dnStrSpace();
    const double mstmCRVShieldHalfLength     =  pSTMShieldPipeParams.dnStrWallHalflength();
    const double mstmCRVShieldHalfWidth      =  pSTMMagnetParams.xHalfLength();
    const double mstmCRVShieldHalfHeight     =  pSTMMagnetParams.yHalfLength();


    G4ThreeVector mstmCRVShieldPositionInMu2e   = stmMagnetPositionInMu2e + G4ThreeVector(0.0,0.0, -pSTMMagnetParams.zHalfLength()-mstmCRVShieldDnStrSpace-mstmCRVShieldHalfLength);
    G4ThreeVector mstmCRVShieldPositionInParent = mstmCRVShieldPositionInMu2e - parentCenterInMu2e;

    // Make the box for the collimator wall
    G4Box* crvShieldBox = new G4Box("crvShieldBox",mstmCRVShieldHalfWidth,mstmCRVShieldHalfHeight,mstmCRVShieldHalfLength);
    //Make the tube for the hole
    G4Tubs *crvShieldHole = new G4Tubs( "crvShieldHole", 0.0, pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(), mstmCRVShieldHalfLength+1.0, 0.0, CLHEP::twopi );

    // create wall with hole
    VolumeInfo crvshield;
    crvshield.name = "crvshield";
    crvshield.solid = new G4SubtractionSolid(crvshield.name,crvShieldBox,crvShieldHole,0,G4ThreeVector(0.0,0.0,0.0));

    //Make the tube to shield CRV
    const double crvShieldTubeHalfLength = pSTMShieldPipeParams.pipeHalfLength();
    G4Tubs *crvShieldTubeTemp = new G4Tubs( "crvShieldTubeTemp", pSTMShieldPipeParams.radiusIn(), pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(), crvShieldTubeHalfLength+mstmCRVShieldHalfLength-0.5*pipeCrack, 0.0, CLHEP::twopi );

    G4ThreeVector mstmCRVShieldTubePositionInMu2e      = mstmCRVShieldPositionInMu2e + G4ThreeVector(0.0,0.0,-mstmCRVShieldHalfLength-crvShieldTubeHalfLength-0.1);
    G4ThreeVector mstmCRVShieldTubeInPositionInParent  = mstmCRVShieldPositionInParent + G4ThreeVector(0.0,0.0,-crvShieldTubeHalfLength+0.5*pipeCrack-0.1);
    G4ThreeVector mstmCRVShieldTubePositionInParent    = mstmCRVShieldPositionInParent + G4ThreeVector(0.0,0.0,-mstmCRVShieldHalfLength-crvShieldTubeHalfLength+0.5*pipeCrack-0.1);

    CLHEP::Hep3Vector vdSTM_UpStrPositionWRTcrvShieldTube           = vdSTM_UpStrPositionInMu2e           - mstmCRVShieldTubePositionInMu2e;
    CLHEP::Hep3Vector vdDSNeutronShieldExitPositionWRTcrvShieldTube = vdDSNeutronShieldExitPositionInMu2e - mstmCRVShieldTubePositionInMu2e;

    G4SubtractionSolid *crvShieldTubeTemp2 = new G4SubtractionSolid("crvShieldTubeTemp2",crvShieldTubeTemp, aDiskVDDSNeutronShieldExitTub, 0, vdDSNeutronShieldExitPositionWRTcrvShieldTube);

    VolumeInfo crvlinershieldtube;
    crvlinershieldtube.name = "crvlinershieldtube";
    crvlinershieldtube.solid = new G4SubtractionSolid(crvlinershieldtube.name,crvShieldTubeTemp2, aDiskVDSTM_UpStrTub, 0, vdSTM_UpStrPositionWRTcrvShieldTube);

    VolumeInfo crvsteelshieldtube;
    crvsteelshieldtube.name = "crvsteelshieldtube";
    crvsteelshieldtube.solid = new G4Tubs(crvsteelshieldtube.name , pSTMShieldPipeParams.radiusIn()+pSTMShieldPipeParams.linerWidth(),  pSTMShieldPipeParams.radiusOut(), crvShieldTubeHalfLength-0.5*pipeCrack, 0.0, CLHEP::twopi );

    if (pSTMShieldPipeParams.build()){
      finishNesting(crvshield,
                    findMaterialOrThrow(pSTMShieldPipeParams.material()),
                    0,
                    mstmCRVShieldPositionInParent,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
      finishNesting(crvlinershieldtube,
                    findMaterialOrThrow(pSTMShieldPipeParams.materialLiner()),
                    0,
                    mstmCRVShieldTubeInPositionInParent,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
      finishNesting(crvsteelshieldtube,
                    findMaterialOrThrow(pSTMShieldPipeParams.material()),
                    0,
                    mstmCRVShieldTubePositionInParent,
                    parentInfo.logical,
                    0,
                    STMisVisible,
                    G4Colour::Magenta(),
                    STMisSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck);
    }



		//wyq shielding for the Ge detector 20170814
		//wyq
		if(pSTMShieldingParams.build()){
			const double stmShieldingHalfLength = pSTMShieldingParams.halfLength();
			const double stmShieldingHalfWidth  = pSTMShieldingParams.halfWidth();
			const double stmShieldingHalfHeight = pSTMShieldingParams.halfHeight();
			const double stmShieldingHoleHalfLength = pSTMShieldingParams.holeHalfLength()+0.001;//break surface
			const double stmShieldingHoleHalfWidth  = pSTMShieldingParams.holeHalfWidth();
			const double stmShieldingHoleHalfHeight = pSTMShieldingParams.holeHalfHeight();
			const double stmShieldingHoleZoffset = pSTMShieldingParams.holeZoffset();
			const double stmShieldingHoleYoffset = pSTMShieldingParams.holeYoffset();
			const double stmShieldingHoleXoffset = pSTMShieldingParams.holeXoffset();
			G4Box* boxShielding1=new G4Box("boxShielding1",stmShieldingHalfWidth,stmShieldingHalfHeight,stmShieldingHalfLength);//wyq
			G4Box* boxShieldingHole=new G4Box("boxShieldingHole",stmShieldingHoleHalfWidth,stmShieldingHoleHalfHeight,stmShieldingHoleHalfLength);//wyq
			G4SubtractionSolid *boxShielding=new G4SubtractionSolid("boxShielding",boxShielding1,boxShieldingHole,0,G4ThreeVector(stmShieldingHoleXoffset,stmShieldingHoleYoffset,stmShieldingHoleZoffset));       //wyq
			//G4SubtractionSolid *boxShielding=new G4SubtractionSolid("boxShielding",boxShielding1,boxShieldingHole,0,G4ThreeVector(0.0,0.0,stmShieldingHoleZshift));       //wyq
			G4ThreeVector stmShieldingPositionInMu2e1=pSTMShieldingParams.originInMu2e();//G4ThreeVector(-3904,0,40900-stmShieldingHalfLength);//wyq
			G4ThreeVector stmShieldingPositionInParent1 = stmShieldingPositionInMu2e1 - parentCenterInMu2e;//wyq
			VolumeInfo shielding;//wyq//wyq
			shielding.name = "shielding";//wyq//wyq
			shielding.solid=boxShielding; //wyq//wyq
			finishNesting(shielding,//wyq//wyq
					findMaterialOrThrow(pSTMShieldingParams.material()),//wyq//wyq
					0,//wyq//wyq
					stmShieldingPositionInParent1,//wyq//wyq
					parentInfo.logical,//wyq//wyq
					0,//wyq//wyq
					true,//wyq//wyq
					G4Colour::Magenta(),//wyq//wyq
					false,//wyq//wyq
					forceAuxEdgeVisible,//wyq//wyq
					placePV,//wyq//wyq
					doSurfaceCheck);//wyq//wyq
		} // end of constructSTM;
	}
}
