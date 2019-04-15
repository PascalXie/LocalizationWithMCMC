int alice()
{
	TH2D * hMCMC = new TH2D("","",100,0.8,1.2,100,1.8,2.2);
	//
	// step 1 : reading data
	//
	ifstream read("data_samples.txt");

	if(read.fail())
	{
		cout<<"File dose not exist."<<endl;
		return 0;
	}

    double x = 0;
    double y = 0;
    double I = 0;

	while(!read.eof())
    {
		read>>x>>y>>I;

		if(read.eof()) break;

        cout<<x<<" "<<y<<" "<<I<<endl;

        hMCMC->Fill(x,y);

    }

	read.close();

	TCanvas * c1 = new TCanvas("c1");
	c1->cd();
	hMCMC->Draw("colz");
	
	hMCMC->GetXaxis()->SetTitle("x / m");
	hMCMC->GetXaxis()->SetLabelFont(12);
	hMCMC->GetXaxis()->SetLabelSize(0.05);
	hMCMC->GetXaxis()->SetTitleSize(0.05);
	hMCMC->GetXaxis()->SetTitleFont(22);
	hMCMC->GetXaxis()->SetTitleOffset(0.9);
	hMCMC->GetYaxis()->SetTitle("y / m");
	hMCMC->GetYaxis()->SetLabelFont(12);
	hMCMC->GetYaxis()->SetLabelSize(0.05);
	hMCMC->GetYaxis()->SetTitleSize(0.05);
	hMCMC->GetYaxis()->SetTitleOffset(0.9);
	hMCMC->GetYaxis()->SetTitleFont(22);

	hMCMC->SetStats(kFALSE);

    return 0;
}
