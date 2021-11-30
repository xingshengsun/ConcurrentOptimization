$------------------------------------------------+
*KEYWORD
$------------------------------------------------+
$## Note:
$    Being run from subdirectory and all paths
$    need to be absolute to THAT location
$------------------------------------------------+
*INCLUDE
$------------------------------------------------+
../../BallElement_D0p747.k 
$------------------------------------------------+
*INCLUDE
$------------------------------------------------+
PlateElement.k
$------------------------------------------------+
*MAT_RIGID
$------------------------------------------------+
$ Steel Projectile
$MID RO E PR N COUPLE M ALIAS/RE
1,7.83,210e10,0.30,,,,
$CMO CON1 CON2
,,
$LCO/A1 A2 A3 V1 V2 V3
,,,,,
$------------------------------------------------+
*MAT_JOHNSON_COOK
$------------------------------------------------+
$ AZ31B plate
$MID,rho,G,E,PR,DTF,VP,Rateop
2,1.77,1.70e11,4.5e11,0.35,,1,0
$  A  B  N C M TM  TR  EPSO
&A,1773984e3,0.294247,0.012701,1.537086,905,298,1e-3
$  CP  PC  SPALL IT D1 D2  D3  D4
1.75e7,-1.5e10,2.0,,-0.381831,&D2,-0.492246,0.479636
$  D5  C2/P/XNP EROD EFMIN
7.754436,,0.0,
$------------------------------------------------+
$ Note: The specific heat for AZ31B is 1.05e7
$ The conversion of plastic work to heat is 0.6
$ So the CP value is set to 1.05e7/0.6=1.75e7
$------------------------------------------------+
*EOS_GRUNEISEN
$------------------------------------------------+
$ AZ31B plate
$EOSID C S1 S2 S3 GAMAO A E0
2,4.52E5,1.242,,,1.54,1.54,0
$V0
1
$------------------------------------------------+
*MAT_ADD_EROSION
$------------------------------------------------+
$ AZ31B plate
$MID EXCL MXPRES MNEPS EFFEPS VOLEPS NUMFIP NCS
2,,,,,,,
$MNPRES SIGP1 SIGVM MXEPS EPSSH SIGTH IMPULSE FAILTM
,,,0.5,,,,
$------------------------------------------------+
*MAT_HYPERELASTIC_RUBBER
$------------------------------------------------+
$ Polyurea
$MID RO PR N NV G SIGF REF
3,1.07,0.485,0,4,0.0224e10,,
$C10 C01 C11 C20 C02 C30
&C10,0.7e7,-0.03e7,-0.02e7,0.001e7,0.002e7
$GI BETAI KI BETAKI
0.0189e10,2.158,,
$GI BETAI KI BETAKI
0.0378e10,1.5608e4,,
$GI BETAI KI BETAKI
0.0805e10,8.5985e6,,
$GI BETAI KI BETAKI
0.0973e10,1.3659e9,,
$------------------------------------------------+
*MAT_ADD_EROSION
$------------------------------------------------+
$ Polyurea
$MID EXCL MXPRES MNEPS EFFEPS VOLEPS NUMFIP NCS
3,,,,,,,
$MNPRES SIGP1 SIGVM MXEPS EPSSH SIGTH IMPULSE FAILTM
,&sig_max,,1.5,,,,
$------------------------------------------------+
*PART
$------------------------------------------------+
Projectile
$PID,SECID,MID,EOSID,HGID,GRAV,ADPOPT,TMID
1000000,1,1,,,,,
$------------------------------------------------+
*PART
$------------------------------------------------+
Plate_top_periphery
$PID,SECID,MID,EOSID,HGID,GRAV,ADPOPT,TMID
3000001,1,2,2,1,,,
$------------------------------------------------+
*PART
$------------------------------------------------+
Plate_top_periphery
$PID,SECID,MID,EOSID,HGID,GRAV,ADPOPT,TMID
3000002,1,2,2,1,,,
$------------------------------------------------+
*PART
$------------------------------------------------+
Plate_top_center
$PID,SECID,MID,EOSID,HGID,GRAV,ADPOPT,TMID
3000003,1,2,2,1,,,
$------------------------------------------------+
*PART
$------------------------------------------------+
Plate_bottom_periphery
$PID,SECID,MID,EOSID,HGID,GRAV,ADPOPT,TMID
3000004,1,3,,2,,,
$------------------------------------------------+
*PART
$------------------------------------------------+
Plate_bottom_periphery
$PID,SECID,MID,EOSID,HGID,GRAV,ADPOPT,TMID
3000005,1,3,,2,,,
$------------------------------------------------+
*PART
$------------------------------------------------+
Plate_bottom_center
$PID,SECID,MID,EOSID,HGID,GRAV,ADPOPT,TMID
3000006,1,3,,2,,,
$------------------------------------------------+
*SECTION_SOLID
$------------------------------------------------+
$SECID,ELFORM,AET
1,1,
$------------------------------------------------+
*CONTACT_ERODING_SURFACE_TO_SURFACE
$------------------------------------------------+
$ Between projectile and top center plate
$SSID MSID SSTYP MSTYP SBOXID MBOXID SPR MPR
3000003,1000000,3,3,,,,
$FS FD DC VC VDC PENCHK BT DT
0.0,0.0,,,,,,
$SFS SFM SST MST SFST SFMT FSF VSF
1.0,1.0,,,,,,
$ISYM EROSOP IADJ
,1,1
$------------------------------------------------+
*CONTACT_ERODING_SURFACE_TO_SURFACE
$------------------------------------------------+
$ Between projectile and bottom center plate
$SSID MSID SSTYP MSTYP SBOXID MBOXID SPR MPR
3000006,1000000,3,3,,,,
$FS FD DC VC VDC PENCHK BT DT
0.0,0.0,,,,,,
$SFS SFM SST MST SFST SFMT FSF VSF
10.0,10.0,,,,,,
$ISYM EROSOP IADJ
,1,1
$------------------------------------------------+
*BOUNDARY_SPC_SET
$------------------------------------------------+
$Symmetry condition for plate in middle XZ plane
$setID,CID,x,y,z,rx,ry,rz
2,0,0,1,0,1,0,1
$------------------------------------------------+
*BOUNDARY_SPC_SET
$------------------------------------------------+
$Symmetry condition for plate in middle YZ plane
$setID,CID,x,y,z,rx,ry,rz
1,0,1,0,0,0,1,1
$------------------------------------------------+
*BOUNDARY_SPC_SET
$------------------------------------------------+
$Symmetry condition for ball in middle XZ plane
$setID,CID,x,y,z,rx,ry,rz
1000000,0,0,1,0,1,0,1
$------------------------------------------------+
*BOUNDARY_SPC_SET
$------------------------------------------------+
$Symmetry condition for ball in middle YZ plane
$setID,CID,x,y,z,rx,ry,rz
1000001,0,1,0,0,0,1,1
$------------------------------------------------+
*INITIAL_VELOCITY_GENERATION
$------------------------------------------------+
$ID STYP OMEGA VX VY VZ IVATN ICID
1000000,2,0,0,0,&vi
$XC YC ZC NX NY NZ PHASE IRIGID
,,,,,,,
$------------------------------------------------+
*INITIAL_TEMPERATURE_SET
$------------------------------------------------+
0,298.15
$------------------------------------------------+
*HOURGLASS
$------------------------------------------------+
$HGID IHQ QM IBQ Q1 Q2 QB/VDC QW
1,6,,,,,,
$------------------------------------------------+
*HOURGLASS
$------------------------------------------------+
$HGID IHQ QM IBQ Q1 Q2 QB/VDC QW
2,4,,,,,,
$------------------------------------------------+
*CONTROL_SOLUTION
$------------------------------------------------+
0
$------------------------------------------------+
*CONTROL_TIMESTEP
$------------------------------------------------+
$DTINIT TSSFAC ISDO TSLIMT DT2MS LCTM ERODE MS1ST
1e-12,0.5,,,,,,
$------------------------------------------------+
*CONTROL_TERMINATION
$------------------------------------------------+
&t_max,20000
$------------------------------------------------+
*TERMINATION_BODY
$------------------------------------------------+
$PID STOP MAXC MINC
1000000,3,,-3.0
$------------------------------------------------+
*CONTROL_ENERGY
$------------------------------------------------+
$HGEN RWEN SLNTEN RYLEN
2,2,2,2
$------------------------------------------------+
*DATABASE_BINARY_D3PLOT
$------------------------------------------------+
1e-6
$------------------------------------------------+
*DATABASE_EXTENT_BINARY
$------------------------------------------------+
$NEIPH NEIPS MAXINT STRFLG SIGFLG EPSFLG RLTFLG ENGFLG
6,,,1,,,,
$CMPFLG IEVERP BEAMIP DCOMP SHGE STSSZ N3THDT IALEMAT
,,,,,,,
$NINTSLD PKP_SEN SCLP HYDRO MSSCL THERM INTOUT NODOUT
,,,,,,,
$DTDT RESPLT NEIPB
,,
$------------------------------------------------+
*DATABASE_MATSUM
$------------------------------------------------+
$DT BINARY LCUR IOOPT OPTION1 OPTION2 OPTION3 OPTION4
1e-7,1,,,,,,
$------------------------------------------------+
*DATABASE_GLSTAT
$------------------------------------------------+
$DT BINARY LCUR IOOPT OPTION1 OPTION2 OPTION3 OPTION4
1e-7,1,,,,,,
$------------------------------------------------+
*DATABASE_HISTORY_NODE_ID
$------------------------------------------------+
$ID HEADING
1000003, BallCenter
$------------------------------------------------+
*DATABASE_NODOUT
$------------------------------------------------+
$DT BINARY LCUR IOOPT OPTION1 OPTION2 OPTION3 OPTION4
1e-7,1
$------------------------------------------------+
*End
$------------------------------------------------+
