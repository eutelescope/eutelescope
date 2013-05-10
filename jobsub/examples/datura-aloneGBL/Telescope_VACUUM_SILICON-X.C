void Telescope_VACUUM_SILICON-X() {
   gSystem->Load("libGeom");
   new TGeoManager("Telescope", "v0.1");

   Double_t dx,dy,dz;
   Double_t dx1, dx2, dy1, dy2;
   Double_t vert[20], par[20];
   Double_t theta, phi, h1, bl1, tl1, alpha1, h2, bl2, tl2, alpha2;
   Double_t twist;
   Double_t origin[3];
   Double_t rmin, rmax, rmin1, rmax1, rmin2, rmax2;
   Double_t r, rlo, rhi;
   Double_t phi1, phi2;
   Double_t a,b;
   Double_t point[3], norm[3];
   Double_t rin, stin, rout, stout;
   Double_t thx, phx, thy, phy, thz, phz;
   Double_t alpha, theta1, theta2, phi1, phi2, dphi;
   Double_t tr[3], rot[9];
   Double_t z, density, radl, absl, w;
   Double_t lx,ly,lz,tx,ty,tz;
   Double_t xvert[50], yvert[50];
   Double_t zsect,x0,y0,scale0;
   Int_t nel, numed, nz, nedges, nvert;
   TGeoBoolNode *pBoolNode = 0;

   // MATERIALS, MIXTURES AND TRACKING MEDIA
// Material: material_World_VACUUM
   a       = 0.000000;
   z       = 0.000000;
   density = 0.000000;
   radl    = 1000000000000000000000000000000.000000;
   absl    = 1000000000000000000000000000000.000000;
   pMat1 = new TGeoMaterial("material_World_VACUUM", a,z,density,radl,absl);
   pMat1->SetIndex(0);
// Medium: medium_World_VACUUM
   numed   = 0;  // medium number
   par[0]  = 0.000000; // isvol
   par[1]  = 0.000000; // ifield
   par[2]  = 0.000000; // fieldm
   par[3]  = 0.000000; // tmaxfd
   par[4]  = 0.000000; // stemax
   par[5]  = 0.000000; // deemax
   par[6]  = 0.000000; // epsil
   par[7]  = 0.000000; // stmin
   pMed1 = new TGeoMedium("medium_World_VACUUM", numed,pMat1, par);
// Material: material_Sensor_SILICON
   a       = 28.085500;
   z       = 14.000000;
   density = 2.330000;
   radl    = 9.349607;
   absl    = 45.753206;
   pMat2 = new TGeoMaterial("material_Sensor_SILICON", a,z,density,radl,absl);
   pMat2->SetIndex(1);
// Medium: medium_Sensor_SILICON
   numed   = 0;  // medium number
   par[0]  = 0.000000; // isvol
   par[1]  = 0.000000; // ifield
   par[2]  = 0.000000; // fieldm
   par[3]  = 0.000000; // tmaxfd
   par[4]  = 0.000000; // stemax
   par[5]  = 0.000000; // deemax
   par[6]  = 0.000000; // epsil
   par[7]  = 0.000000; // stmin
   pMed2 = new TGeoMedium("medium_Sensor_SILICON", numed,pMat2, par);

   // TRANSFORMATION MATRICES
   // Translation: matrix_Sensor1
   dx = 0.000000;
   dy = 0.000000;
   dz = 0.000000;
   TGeoTranslation *pMatrix2 = new TGeoTranslation("matrix_Sensor1",dx,dy,dz);
   // Translation: matrix_Sensor2
   dx = 150.000000;
   dy = 0.000000;
   dz = 0.000000;
   TGeoTranslation *pMatrix4 = new TGeoTranslation("matrix_Sensor2",dx,dy,dz);
   // Translation: matrix_Sensor3
   dx = 300.000000;
   dy = 0.000000;
   dz = 0.000000;
   TGeoTranslation *pMatrix5 = new TGeoTranslation("matrix_Sensor3",dx,dy,dz);
   // Translation: matrix_Sensor4
   dx = 450.000000;
   dy = 0.000000;
   dz = 0.000000;
   TGeoTranslation *pMatrix6 = new TGeoTranslation("matrix_Sensor4",dx,dy,dz);
   // Translation: matrix_Sensor5
   dx = 600.000000;
   dy = 0.000000;
   dz = 0.000000;
   TGeoTranslation *pMatrix7 = new TGeoTranslation("matrix_Sensor5",dx,dy,dz);
   // Translation: matrix_Sensor6
   dx = 750.000000;
   dy = 0.000000;
   dz = 0.000000;
   TGeoTranslation *pMatrix8 = new TGeoTranslation("matrix_Sensor6",dx,dy,dz);
   // Shape: Box_World type: TGeoBBox
   dx = 5000.000000;
   dy = 5000.000000;
   dz = 5000.000000;
   TGeoShape *pBox_World_1 = new TGeoBBox("Box_World", dx,dy,dz);
   // Volume: volume_World
   pvolume_World_3f27e48 = new TGeoVolume("volume_World",pBox_World_1, pMed1);
   pvolume_World_3f27e48->SetLineColor(4);
   pvolume_World_3f27e48->SetLineWidth(3);
   pvolume_World_3f27e48->SetVisLeaves(kTRUE);

   // SET TOP VOLUME OF GEOMETRY
   gGeoManager->SetTopVolume(pvolume_World_3f27e48);

   // SHAPES, VOLUMES AND GEOMETRICAL HIERARCHY
   // Shape: Box_Sensor type: TGeoBBox
   dx = 0.025000;
   dy = 10.750000;
   dz = 6.850000;
   TGeoShape *pBox_Sensor_2 = new TGeoBBox("Box_Sensor", dx,dy,dz);
   // Volume: volume_Sensor1
   pvolume_Sensor1_4202828 = new TGeoVolume("volume_Sensor1",pBox_Sensor_2, pMed2);
   pvolume_Sensor1_4202828->SetVisLeaves(kTRUE);
   pvolume_World_3f27e48->AddNode(pvolume_Sensor1_4202828, 1, pMatrix2);
   // Volume: volume_Sensor2
   pvolume_Sensor2_4202718 = new TGeoVolume("volume_Sensor2",pBox_Sensor_2, pMed2);
   pvolume_Sensor2_4202718->SetVisLeaves(kTRUE);
   pvolume_World_3f27e48->AddNode(pvolume_Sensor2_4202718, 2, pMatrix4);
   // Volume: volume_Sensor3
   pvolume_Sensor3_42028b0 = new TGeoVolume("volume_Sensor3",pBox_Sensor_2, pMed2);
   pvolume_Sensor3_42028b0->SetVisLeaves(kTRUE);
   pvolume_World_3f27e48->AddNode(pvolume_Sensor3_42028b0, 3, pMatrix5);
   // Volume: volume_Sensor4
   pvolume_Sensor4_4202938 = new TGeoVolume("volume_Sensor4",pBox_Sensor_2, pMed2);
   pvolume_Sensor4_4202938->SetVisLeaves(kTRUE);
   pvolume_World_3f27e48->AddNode(pvolume_Sensor4_4202938, 4, pMatrix6);
   // Volume: volume_Sensor5
   pvolume_Sensor5_42029c0 = new TGeoVolume("volume_Sensor5",pBox_Sensor_2, pMed2);
   pvolume_Sensor5_42029c0->SetVisLeaves(kTRUE);
   pvolume_World_3f27e48->AddNode(pvolume_Sensor5_42029c0, 5, pMatrix7);
   // Volume: volume_Sensor6
   pvolume_Sensor6_4202a48 = new TGeoVolume("volume_Sensor6",pBox_Sensor_2, pMed2);
   pvolume_Sensor6_4202a48->SetVisLeaves(kTRUE);
   pvolume_World_3f27e48->AddNode(pvolume_Sensor6_4202a48, 6, pMatrix8);

   // CLOSE GEOMETRY
   gGeoManager->CloseGeometry();
   gGeoManager->GetTopVolume()->Draw();
}
