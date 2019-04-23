/***************************************************************************
 *  How to use ProMC files from HepSim, and how to build anti-KT jets 
 *  S.Chekanov (ANL) chekanov@anl.gov
 *  A library for HEP events storage and processing based on Google's PB   
 *  The project web site: http://atlaswww.hep.anl.gov/hepsim/
****************************************************************************/

#include<iostream>
#include<fstream>
#include<stdlib.h>

// check directory
#include <sys/stat.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH1D.h>
#include "TMath.h"
#include"time.h"
#include <dirent.h>
#include <string>
#include <vector>
#include <map>
struct stat sb;

// fastjet
#include "fjcore.hh"
// promc
#include "ProMC.pb.h"
#include "ProMCBook.h"

#define pi  3.141592654
#define pi2 6.283185308

using namespace std;
using namespace promc;

void kcluster (int nclusters, int nrows, int ncolumns,
               double**, int**, double weight[], int,
               int, char, char,
               int clusterid[], double*, int*);

double clusterdistance (int nrows, int ncolumns, double** data,
                        int** mask, double weight[], int n1, int n2, int index1[], int index2[],
                        char dist, char method, int transpose);

void getclustermean(int nclusters, int nrows, int ncolumns,
                    double** data, int** mask, int clusterid[], double** cdata, int** cmask,
                    int transpose);


// find all files inside a directory
std::vector<std::string> open(std::string path = ".") {
	std::vector<std::string> files;
	if (stat(path.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)) {
		DIR*    dir;
		dirent* pdir;
		dir = opendir(path.c_str());
		while (pdir = readdir(dir)) {
			string rfile=path + "/" + pdir->d_name;
			if (rfile.find(".promc") != std::string::npos)
				files.push_back(rfile);
		}
	} else {
		cout << "Directory " << path << " does not exist!!" << endl;
		std::exit(0);
	}
	return files;
}


// main example
int main(int argc, char **argv)
{
	// total events
	int ntot=0; // total events
	int nfiles=0; // total files


	string dir="./data/";
	std::vector<std::string> files = open(dir);

	double cross=0;
	double xcross=0;
	int nselect=0;

	string outputfile="output.root";
	cout << "\n -> Output file is =" << outputfile << endl;
	TFile * RootFile = new TFile(outputfile.c_str(), "RECREATE", "Histogram file");
	//  RootFile->SetCompressionLevel(0);

	TH1D * h_debug = new TH1D("debug", "debug", 10, 0, 10.);
	TH1D * h_cross = new TH1D("cross", "cross,events,lumi", 5, 0, 5.);
	TH1D * h_pt = new TH1D("jet_pt", "pt",100,0,1000);
	TH1D * h_eta = new TH1D("jet_eta", "eta", 40, -10, 10);

	// fastjet
	const double R = 0.4;
        const int nclusters = 6;
        const int ncols = 2;
        int ifound = 0;
        double error;
        const int transpose = 0;
        const char dist = 'r';
        const char method = 'a';
        const int npass = 300; // number of passes for k-means 
        const double EtaMax=3.0; //  max eta

	fjcore::JetDefinition jet_def(fjcore::antikt_algorithm, R);
	std::map<int,int> chargemap; // pad that keeps charge*3

        cout << "Merge particles into " << nclusters << " jets/clusters" << endl;
        cout << "Jet cone size is     " << R << endl;

	for(unsigned int m=0; m < files.size(); m++){
		string Rfile=files[m];
		ProMCBook*  epbook = new ProMCBook(Rfile.c_str(),"r");

		cout << "\n\n Start to read.." << endl;
		// get the version number
		int  h=epbook->getVersion();
		if (m==0) cout << "Version = " << h << endl;
		// get the description of this file
		string d=epbook->getDescription();
		if (m==0) cout << "Description = " << d << endl;
		int  nev=epbook->getEvents();
		cout << "Events = " << nev  << endl;
		// get the header file with units, cross sections etc.
		ProMCHeader header = epbook->getHeader();


		// this part reads the header information with particle data.
		// you can access names, id, charge*3, and other useful
		// information from PDG table. As an example, we make a map
		// that keeps charge*3.
		if (m==0) {
			for (int jj=0; jj<header.particledata_size(); jj++){
				ProMCHeader_ParticleData data= header.particledata(jj);
				int charge3=data.charge();
				int id=data.id();
				double mass=data.mass();
				string name=data.name();
				cout << "name=" << name << " mass=" << mass << " charge=" << charge3 << endl;
				chargemap[id]=charge3;
			}
		}


		// here are the units
		double kEV=(double)(header.momentumunit());
		//double kLe=(double)(header.lengthunit());


		// loop over all events
		for (int nn=0; nn<nev; nn++){
			if (epbook->next() !=0) continue;
			ProMCEvent eve = epbook->get();

			// get truth information
			ProMCEvent_Particles  *pa=eve.mutable_particles();
			h_debug->Fill(1.0);
			ntot++;

			double xlumi=0;
			if (m>0) {
				xlumi=(double)ntot/cross; // lumi so far
				h_debug->Fill(2.0);
			}

			if (ntot%1000==0)
				cout <<  " # Events=" << ntot  << " X-cross="<< cross << " pb" << " Lumi so far=" << xlumi/1000.0 << " fb-1" << endl;

			vector<fjcore::PseudoJet> avec;



			// fill stable and no neutrino
			for (int j=0; j<pa->pdg_id_size(); j++){
				if (pa->status(j)!=1) continue;
				if (abs(pa->pdg_id(j)==12) || abs(pa->pdg_id(j))==14 || abs(pa->pdg_id(j))==16 ) continue;
				double px= pa->px(j)/kEV;
				double py= pa->py(j)/kEV;
				double pz= pa->pz(j)/kEV;
				double e= pa->energy(j)/kEV;
				double pt=sqrt(px*px+py*py);
				double eta=-log(tan(atan2(pt,(double)pz)/2));
				// int charge=chargemap[pa->pdg_id(j)]; // get charge

				if ( pt < 0.3)                   continue;
				if ( fabs(eta)> EtaMax )         continue;
				avec.push_back( fjcore::PseudoJet(px,py,pz,e) );
			}


			fjcore::ClusterSequence clust_seq(avec, jet_def);
                        vector<fjcore::PseudoJet> exclusive_jets = clust_seq.exclusive_jets_up_to( nclusters );
			// sort in PT
			vector<fjcore::PseudoJet> sorted_jets = sorted_by_pt(exclusive_jets);

			for (unsigned int k = 0; k<sorted_jets.size(); k++) {
				double eta=sorted_jets[k].pseudorapidity();
				double phi=sorted_jets[k].phi();
				double pt = sorted_jets[k].perp();
				double e = sorted_jets[k].e();
				h_pt->Fill(pt);
				h_eta->Fill(eta);
			} // end loop


			// now k-means clustering
			int nrows = avec.size();
			double** data = (double**)malloc(nrows*sizeof(double*) );
			int** mask = (int**)malloc(nrows*sizeof(int*));
			if(data == NULL || mask == NULL )
			{
				printf("Not enough memory\n");
				exit(1);
			}
			for (int i=0; i<nrows; i++)
			{ data[i] = (double*)malloc(ncols*sizeof(double));
				mask[i] = (int*)malloc(ncols*sizeof(int));
			}


			for  (int i=0; i<nrows; i++ ) {
				data[i][0]=0.0;
				data[i][1]=0.0;

				fjcore::PseudoJet pp = avec[i];
				if (pp.perp()>0) {
					data[i][0]=pp.pseudorapidity();
					data[i][1]=pp.phi();
					if (data[i][1]<0) data[i][1]=data[i][1]+pi2;
				}

				// // it is important that each metric contribute equally to the total distance.
				// // shift to positive and normilize
				data[i][0]=(EtaMax+data[i][0])/(2*EtaMax);
				data[i][1]=data[i][1]/pi2;

				for  (int j=0; j<ncols; j++ )  {
					mask[i][j]=1;
					if (data[i][j] == 0.0) mask[i][j]=0;
				}
			}


                        // add energy weights
			double* weight = (double*)malloc(ncols*sizeof(double));
			for (int i= 0; i<ncols; i++) {
                                               fjcore::PseudoJet pp = avec[i];
                                               weight[i] = 1.0; // pp.e();

                        }
			int* clusterid = (int*)malloc(nrows*sizeof(int));

			// do k-means
			// printf("%d  pass of the K-means algorithm (result should not change)\n",npass);
			kcluster(nclusters,nrows,ncols,data,mask,weight,transpose,npass,method,dist,
			         clusterid, &error, &ifound);
			//   printf ("k-means clustering:  Solution found %d times;\n ", ifound);
			//   printf ("k-means clustering:  Solution corresponds to  %12.5lf sum of distances;\n", error);


			double distance;
			int** index;
			int* count;
			int i;

			/*
			  printf ("Cluster assignments:\n");
			  for (i = 0; i < nrows; i++)
			    printf ("Point %d: cluster %d\n", i, clusterid[i]);
			  printf ("\n");
			*/

			//  printf ("------- Distance between clusters:\n");
			index = (int**)malloc(nclusters*sizeof(int*));
			count = (int*)malloc(nclusters*sizeof(int));
			for (i = 0; i < nclusters; i++) count[i] = 0;
			for (i = 0; i < nrows; i++) count[clusterid[i]]++;
			for (i = 0; i < nclusters; i++) index[i] = (int*)malloc(count[i]*sizeof(int));
			for (i = 0; i < nclusters; i++) count[i] = 0;
			for (i = 0; i < nrows; i++)
			{ int id = clusterid[i];
				index[id][count[id]] = i;
				count[id]++;
			}
			distance =
			    clusterdistance(nrows, ncols, data, mask, weight, count[0], count[1],
			                    index[0], index[1], dist, method, 0);
			//  printf("Distance between 0 and 1: %7.3f\n", distance);


			//free memory
			for (int j= 0; j<nclusters; j++) free(index[j]);
			free(index);
			free(count);


			double** cdata = (double**)malloc(nclusters*sizeof(double*));
			int** cmask = (int**)malloc(nclusters*sizeof(int*));
			for (i = 0; i < nclusters; i++)
			{ cdata[i] = (double*)malloc(ncols*sizeof(double));
				cmask[i] = (int*)malloc(ncols*sizeof(int));
			}



                        if (ntot%10==0) {
                        cout <<  " # Events=" << ntot  << endl;
			printf ("\n");
			printf ("------- k-means cluster centroids:\n");
			getclustermean(nclusters, nrows, ncols, data, mask, clusterid,
			               cdata, cmask, 0);
			for (i = 0; i < nclusters; i++) { 


                         printf("k-means cluster     %2d:", i);

                         // back to normal metric 
                         cdata[i][0] = cdata[i][0]*2*EtaMax - EtaMax; 
                         cdata[i][1] = cdata[i][1]*pi2;


                         for (int j = 0; j < ncols; j++) printf("\t%7.3f", cdata[i][j]);
                         printf("\n");
                         }

                         printf ("------- Jet centers:\n"); 
                         for (i = 0; i < nclusters; i++) {
                             printf("jet center position %2d:", i);
                             printf("\t%7.3f", sorted_jets[i].rapidity()); printf("\t%7.3f", sorted_jets[i].phi()); // printf("\t%7.3f", sorted_jets[i].perp());
                             printf("\n");
			}


			printf("\n");


                        }


		    for (int  j= 0; j<nrows; j++) { free(data[j]); free(mask[j]); }
			free(data);
			free(mask);

			free(weight);
			free(clusterid);


			for (int j=0; j<nclusters; j++)
			{ free(cdata[j]);
				free(cmask[j]);
			}
			free(cdata);
			free(cmask);




		} // end event loop


		ProMCStat stat = epbook->getStatistics();
		cross=stat.cross_section_accumulated();
		epbook->close(); // close
		nfiles++;
		xcross=xcross+cross;


	} // end loop over all files



	xcross=xcross/(double)nfiles; // average cross for all files
	cout << "Total events=" << ntot << endl;
	cout << "Total files=" << nfiles << endl;
	cout << "Average cross section for all files = " << xcross << " pb"<< endl;
	double width=h_pt->GetBinWidth(1);
	double lumi=(double)ntot/xcross;
	cout << "Lumi for all files  =" << lumi << " pb-1" << endl;
	double norm=width*lumi;
	h_pt->Scale(1.0/norm);

	h_cross->Fill(0.0,(double)xcross); // 1
	h_cross->Fill(1.0,(double)ntot); // 2
	h_cross->Fill(2.0,(double)lumi); // 3
	h_cross->Fill(3.0,(double)nfiles);
	cout << "Nr of selected jets=" << nselect << endl;

	/*
	// calculate cross sections
	   TAxis *xaxis = h_pt->GetXaxis(); 
	   xaxis->SetTitle("p_{T}(jet) GeV");
	   TAxis *yaxis = h_pt->GetYaxis();
	   yaxis->SetTitle("d #sigma / dp_{T} [pb / GeV]");

	    norm=lumi*(h_pt->GetBinWidth(1));
	    h_pt->Scale(1.0/norm); 
	    TAxis *xaxis1 = h_pt->GetXaxis();
	    xaxis1->SetTitle("p_{T}(jet) GeV");
	    TAxis *yaxis1 = h_pt->GetYaxis();
	    yaxis1->SetTitle("d #sigma / dp_{T} [pb / GeV]");
	*/

	//m_ntuple->Fill();
	RootFile->Write();
	RootFile->Print();
	RootFile->Close();

	cout << "Writing ROOT file "+ outputfile << endl;

	return 0;
}
