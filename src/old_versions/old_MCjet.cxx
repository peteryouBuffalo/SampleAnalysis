// If we have at least two 'c-tagged' jets and two 'b-tagged' jets
    // (which we require for our events), find the jets that match to
    // our gen objects. We do this via dR-matching because there is no
    // variable to match them.
    if (genCjet_list.size() >= 2 && genBjet_list.size() >= 2) {      

      found_MCjets = true;
      
      // Find how the b-jets match the b-quarks.
      std::vector<std::pair<int,float>> jets_idx_dR;
      std::vector<int> chosenIdx;
      // Go through each of the b-partons...
      for (size_t i = 0; i < gen_bs.size(); ++i) {
	
	// ...and check the separation to each possible jet.
        for (size_t j = 0; j < genBjet_list.size(); ++j) {
          float dR = fabs(gen_bs[i].m_lvec.DeltaR(genBjet_list[j].m_lvec));
          //std::cout << "b[" << i << "] | bjet[" << j << "] = " << dR << "\n";
          jets_idx_dR.push_back(std::make_pair(j,dR));
        }

        // Sort the matches via the dR values (second in pair). We sort
	// in ascending order, so our proper choice will be the first option.
        std::sort(jets_idx_dR.begin(), jets_idx_dR.end(), sort_by_second);
        std::pair<int, float> proper_pair = jets_idx_dR[0];

	// Move the proper jet to our list of chosen jets and remove
	// it as an option from the list.
        int idx = proper_pair.first;
        chosenIdx.push_back(idx);
        gen_bjets.push_back(genBjet_list[idx]);
        genBjet_list.erase(genBjet_list.begin() + idx);

        // Fill match information into histograms.
        h_dR_bjets->Fill(proper_pair.second, evtW);
        if (i == 0) h_dR_bjet0->Fill(proper_pair.second, evtW);
        else h_dR_bjet1->Fill(proper_pair.second, evtW);

      }//end-i
      
      // Find how the c-jets match the c-quarks.
      std::vector<std::pair<int,float>> jets_idx_dR2;
      std::vector<int> chosenIdx2;
      // Go through each of the c-partons...
      for (size_t i = 0; i < gen_cs.size(); ++i) {
	
	// ...and check the separation to each possible jet.
        for (size_t j = 0; j < genCjet_list.size(); ++j) {
          float dR = fabs(gen_cs[i].m_lvec.DeltaR(genCjet_list[j].m_lvec));
          //std::cout << "c[" << i << "] | cjet[" << j << "] = " << dR << "\n";
          jets_idx_dR2.push_back(std::make_pair(j,dR));
        }//end-j

        // Sort the matches via the dR values (second in pair). We sort
	// in ascending order, so our proper choice will be the first option.
        std::sort(jets_idx_dR2.begin(), jets_idx_dR2.end(), sort_by_second);
        std::pair<int, float> proper_pair = jets_idx_dR2[0];

	// Move the proper jet to our list of chosen jets and remove
	// it as an option from the list.
        int idx2 = proper_pair.first;
        chosenIdx2.push_back(idx2);
        gen_cjets.push_back(genCjet_list[idx2]);
        genCjet_list.erase(genCjet_list.begin() + idx2);
        
        // Fill match information into histograms.
        h_dR_cjets->Fill(proper_pair.second, evtW);
        if (i == 0) h_dR_cjet0->Fill(proper_pair.second, evtW);
        else h_dR_cjet1->Fill(proper_pair.second, evtW);
      }//end-i
      
      h_evt_VbbHcc->Fill(2.5, evtW);

      // Check to make sure the pT of our Truth jets at least roughly
      // matches the pT of the associated parton.
      bool pT_matches = true;
      for (int i = 0; i < 2; ++i){
        // Check the match of the b-jets
        float pT_jet = gen_bjets[i].Pt();
        float pT_parton = gen_bs[i].Pt();
        float SF = pT_jet/pT_parton;
        h_pT_ratio->Fill(SF, evtW);
        //if (SF <= 0.5) pT_matches = false;

        // Check the match of the c-jets
        pT_jet = gen_cjets[i].Pt();
        pT_parton = gen_cs[i].Pt();
        SF = pT_jet/pT_parton;
        h_pT_ratio->Fill(SF, evtW);
        //if (SF <= 0.5) pT_matches = false; 
      }
      
      if (pT_matches) {

        // From the jets that we've found, reconstruct the proper bosons
        // and fill our desired methods.
        ZObj Z_MCjet(gen_bjets);
        HObj H_MCjet(gen_cjets);

        if (Z_MCjet.M() < 30) {
          for (int i = 0; i < 2; ++i) {
            float dR = jets_idx_dR[chosenIdx[i]].second;
            h_dR_ZUnder30->Fill(dR, evtW); 
          }
        }
        if (H_MCjet.M() < 30) {
          for (int i = 0; i < 2; ++i) {
            float dR = jets_idx_dR2[chosenIdx2[i]].second;
            h_dR_HUnder30->Fill(dR, evtW);
          }
        }
      
        h_VH_MCjet->FillVH(Z_MCjet, H_MCjet, evtW);

      }//end-pT-matches
    }//end-found-jets
    
    /*if (genCjet_list.size() == 1) {
      std::cout << "only ONE cjet found" << std::endl;
      std::cout << ">>> pt = " << genCjet_list[0].Pt() << "\n";
      std::cout << ">>> m  = " << genCjet_list[0].M() << "\n";
    }
    if (genBjet_list.size() == 1) {
      std::cout << "only ONE bjet found" << std::endl;
      std::cout << ">>> pt = " << genBjet_list[0].Pt() << "\n";
      std::cout << ">>> m  = " << genBjet_list[0].M() << "\n";
    }*/
