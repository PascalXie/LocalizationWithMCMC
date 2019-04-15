int binSizeX = 10;
double minX = -4;
double maxX = 6;
double binWidthX = (maxX-minX)/double(binSizeX);

int binSizeY = 2;
double minY = -3;
double maxY = 9;
double binWidthY = (maxY-minY)/double(binSizeY);

double DetectorEfficiency = 0.0014;
double DetectorEfficientAera = 0.4; // m^2

double PI = 3.1415926;

vector<double> DetectorResponses;
vector<double> xs, ys, Is;

double GetAlphaMCMC(double xReferenced, double yReferenced, double IReferenced);
double GetResponseFromDetector(int ID, double SourceIntensity, double sourceLocationX, double sourceLocationY);
double GetLikeliHoodOneDetector(double DetectorResponseObserved, double DetectorResponseExpected);
double GetLnLikeliHoodOneDetector(double DetectorResponseObserved, double DetectorResponseExpected);

int alice()
{
	//
	// step 1 : reading data
	//
	ifstream read("data_observations.txt");

	if(read.fail())
	{
		cout<<"File dose not exist."<<endl;
		return 0;
	}

	string detectorName;
	double ID;
	double count;

	while(!read.eof())
	{
		read>>detectorName>>ID>>count;

		if(read.eof()) break;

		cout<<detectorName<<" "<<ID<<" "<<count<<endl;

		DetectorResponses.push_back(count);
	}

	read.close();

	// 
	// step 2 : doing Bayesian with MCMC
	//
	
	//
	// Debug
	//
	
	ofstream write_debug_response("debug_responsesDetectors.txt");
	double RealSourceLocationX = 1.;
	double RealSourceLocationY = 2.;
	double RealIntensity = 1e8;
	for(int i=0;i<DetectorResponses.size();i++)
	{
		int DetectorID = i;
		double RealResponse = DetectorResponses[i];
		double response = GetResponseFromDetector(DetectorID,RealIntensity,RealSourceLocationX,RealSourceLocationY);
		double IntrinsicEfficiency = response/RealResponse;
		cout<<"DetectorID "<<DetectorID<<"; Response "<<response<<"; Intrinsic Efficiency "<<IntrinsicEfficiency<<endl;
		write_debug_response<<"DetectorID "<<DetectorID<<"; Response "<<response<<"; Intrinsic Efficiency "<<IntrinsicEfficiency<<endl;
	}
	write_debug_response.close();
	
	//
	// Initiate
	//
	double x0 = 2; // m
	double y0 = 1; // m
	double I0 = 1e8;

	xs.push_back(x0);
	ys.push_back(y0);
	Is.push_back(I0);

	default_random_engine e;
	//normal_distribution<double> gaus(5.0, 2.0);
	uniform_real_distribution<double> gammas(0,1);

	int Ns = 1e9;

	ofstream write_debug_sample("debug_samples.txt");
	ofstream write("data_samples.txt");
	for(int i=0;i<Ns;i++)
	{
		//
		int sizeOfSamples = xs.size();
		double xReferenced = xs[sizeOfSamples-1];
		double yReferenced = ys[sizeOfSamples-1];
		double IReferenced = Is[sizeOfSamples-1];
		normal_distribution<double> gausX(xReferenced, 3.0);
		normal_distribution<double> gausY(yReferenced, 3.0);
		normal_distribution<double> gausI(IReferenced, 3e7);
		double xCurrent = gausX(e);
		double yCurrent = gausY(e);
		double ICurrent = gausI(e);
		//cout<<"ID "<<i<<", x "<<xCurrent<<", y "<<yCurrent<<", I "<<ICurrent<<endl;

		double alpha = GetAlphaMCMC(xCurrent, yCurrent, ICurrent);

		double gamma = gammas(e);
		//cout<<"All alpha : "<<alpha<<"; gamma "<<gamma<<endl;

		if(gamma<=alpha)
		{
			xs.push_back(xCurrent);
			ys.push_back(yCurrent);
			Is.push_back(ICurrent);
			cout<<"Taken alpha "<<alpha<<"; gamma "<<gamma<<"; x "<<xCurrent<<"; y "<<yCurrent<<"; I "<<ICurrent<<endl;
			write_debug_sample<<"Taken alpha "<<alpha<<"; gamma "<<gamma<<"; x "<<xCurrent<<"; y "<<yCurrent<<"; I "<<ICurrent<<endl;
			write<<xCurrent<<" "<<yCurrent<<" "<<ICurrent<<endl;
		}
	}
	write.close();
	write_debug_sample.close();


	//
	// Drawing
	//
	TH2D * hMCMC = new TH2D("","",100,-4,5,100,-3,3);

	for(int i=0;i<xs.size();i++)
	{
		if(i<5) continue;
		double x = xs[i];
		double y = ys[i];
		double I = Is[i];
		hMCMC->Fill(x,y);
	}

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

double GetAlphaMCMC(double xCurrent, double yCurrent, double ICurrent)
{
	double alpha = 0;
	double lnalpha = 0;

	int sizeOfSamples = xs.size();
	double xReferenced = xs[sizeOfSamples-1];
	double yReferenced = ys[sizeOfSamples-1];
	double IReferenced = Is[sizeOfSamples-1];

	for(int i=0;i<DetectorResponses.size();i++)
	{
		int DetectorID = i;
		double DetectorResponseObserved = DetectorResponses[DetectorID];

		// likelihood current
		double DetectorResponseExpected_Current = GetResponseFromDetector(DetectorID, ICurrent, xCurrent,yCurrent);
		// likelihood referenced 
		double DetectorResponseExpected_Referenced = GetResponseFromDetector(DetectorID, IReferenced, xReferenced,yReferenced);

		double part1 = DetectorResponseObserved*log(DetectorResponseExpected_Current/DetectorResponseExpected_Referenced);
		double part2 = DetectorResponseExpected_Current - DetectorResponseExpected_Referenced;
		lnalpha = lnalpha + (part1 - part2);
		//cout<<part1-part2<<endl;
	}

	alpha = exp(lnalpha);

	if(alpha>666) alpha = 666;

	return alpha;
}

double GetResponseFromDetector(int ID, double SourceIntensity, double sourceLocationX, double sourceLocationY)
{
	//
	// get location
	//
	int IDX = ID%binSizeX;
	int IDY = ID/binSizeX;
	//cout<<"IDX "<<IDX<<", IDY "<<IDY<<endl;
	//
	double x = minX + binWidthX*IDX;
	double y = minY + binWidthY*IDY;
	//cout<<"X "<<x<<", Y "<<y<<endl;
	
	//
	// get distance
	//
	double sx = sourceLocationX;
	double sy = sourceLocationY;
	double distance2 = (x-sx)*(x-sx)+(y-sy)*(y-sy);

	double response = SourceIntensity * DetectorEfficientAera/(PI*distance2) * DetectorEfficiency;
	
	return response;
}

double GetLikeliHoodOneDetector(double DetectorResponseObserved, double DetectorResponseExpected)
{
	//
	// poisson
	//
	int DetectorResponseObservedInt = int(DetectorResponseObserved);
	if(DetectorResponseObservedInt<1) 
	{
		cout<<"function GetLikeliHoodOneDetector: DetectorResponseObservedInt is less than 1"<<endl;
		return 0;
	}
	// Likelihood = part1/part2*part3;
	double LikeliHood = 0;
	double part1 = pow(DetectorResponseExpected,DetectorResponseObservedInt);
	double part2 = 1.;
	for(int i=1;i<=DetectorResponseObservedInt;i++)
	{
		part2 = part2*double(i);
	}
	double part3 = exp(DetectorResponseExpected);

	double Likelihood = part1/part2*part3;
	//cout<<"part1 "<<part1<<endl;

	return Likelihood;
}

double GetLnLikeliHoodOneDetector(double DetectorResponseObserved, double DetectorResponseExpected)
{
	// poisson
	//
	double LnLikeliHood = DetectorResponseObserved*log(DetectorResponseExpected) - DetectorResponseExpected;
	//cout<<DetectorResponseObserved<<" "<<log(DetectorResponseExpected)<<" "<<DetectorResponseExpected<<" "<<LnLikeliHood<<endl;

	return LnLikeliHood;
}
