void CreateDSChannelMapsFile()
{
	TFile * f = new TFile("DSChannelMaps.root","CREATE");

	GATDataSet ds0(2690); // run numbers in this macro are from ranges in GAT/Apps/DataSetInfo.hh
	MJTChannelMap *chmap0 = ds0.GetChannelMap();
	chmap0->Write("ChannelMapDS0");

	GATDataSet ds1(9639);
	MJTChannelMap *chmap1 = ds1.GetChannelMap();
	chmap1->Write("ChannelMapDS1");

	GATDataSet ds2(15070);
	MJTChannelMap *chmap2 = ds2.GetChannelMap();
	chmap2->Write("ChannelMapDS2");

	GATDataSet ds3(16932);
	MJTChannelMap *chmap3 = ds3.GetChannelMap();
	chmap3->Write("ChannelMapDS3");

	GATDataSet ds4(60001034);
	MJTChannelMap *chmap4 = ds4.GetChannelMap();
	chmap4->Write("ChannelMapDS4");

	GATDataSet ds5(23650);
	MJTChannelMap *chmap5 = ds5.GetChannelMap();
	chmap5->Write("ChannelMapDS5");

  GATDataSet ds6(25800);
  MJTChannelMap *chmap6 = ds6.GetChannelMap();
  chmap6->Write("ChannelMapDS6");

	f->Close();
}
