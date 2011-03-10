/*
	Class which defines an alignment of nucleotide or amino acid sequences.
	Relies on Code, Pair, Sequence and Site classes to implement the algorithm to determine compatibility between sites
*/

import java.util.*;
import java.io.*;
import java.lang.*;

public class Alignment 
{
	private LinkedList<Sequence> sequences; // the sequences making up the alignment
	private LinkedList<Site> informativeSites; // list of the sites found to be informative
	private LinkedList<Site> selectedSites; // list of the sites found to be informative
	private LinkedList<Site> uniqueSites; // list of the informative sites, with ones with redundant character pattern removed
	private LinkedList<Integer> uniquemap; // map logging the number of sites corresponding the the unique character pattern sites
	private LinkedList<Site> allSites; // a list of all the sites
	private LinkedList<Site> binarySites; // a list of the binary informative sites.
	private int size; // number of sequences in the alignment
	private int length; // length of longest sequence in alignment
	private Code code; // Character code ie amino acids or nucleotides
	private RecodeSelector recode_selector; // GUI where user specifies how they want to recode the sequences
	private double ov_comp; // the overall compatibility of the alignment
	private double[] site_comp; // the overall comatibility of each site within the alignment
	public LinkedList<int[][]> mnic_comps;
	private boolean POCcalculated, COMPcalculated;

	/*
		Class Declaration
	*/
	public Alignment(LinkedList<Sequence> sequences, Code c, RecodeSelector r) {
		size = sequences.size();
		if(size == 0)
			length = 0;
		else
			length = Alignment.getShortestSequenceLength(sequences);
		this.sequences = sequences;
		informativeSites = new LinkedList<Site>();
		uniqueSites= new LinkedList<Site>();
		selectedSites= new LinkedList<Site>();
		uniquemap=new LinkedList<Integer>();
		allSites= new LinkedList<Site>();
		binarySites= new LinkedList<Site>();
		code = c;
		recode_selector = r;
		ov_comp=0;
		mnic_comps = new LinkedList<int[][]>();
		site_comp = new double[length];
		POCcalculated = false;
		COMPcalculated = false;
	}
	
	/*
		Retrieval Methods
	*/
	public int informativeTotal(){return informativeSites.size();}
	public int selectedTotal(){return selectedSites.size();}
	public int getSize(){ return size;}
	public int getLength(){return length;}
	public Code getCode(){return code;}
	public RecodeSelector getRecode(){return recode_selector;}
	public LinkedList<Sequence> getSeqs(){return sequences;}
	public double getOvComp(){return ov_comp;}
	public LinkedList<int[][]> getMNIC(){return mnic_comps;}
	public void resetOvComp(){ov_comp=0;}
	public double[] getSiteComp(){return site_comp;}
	public LinkedList<Site> getSites(){	return allSites;}
	
	
	/*
		Returns the length of the shortest string in sequences
	*/
	public static int getShortestSequenceLength(LinkedList<Sequence> sequences) {
		int shortest_length;
		ListIterator<Sequence> iter = sequences.listIterator();
		String tempstring = ((Sequence)iter.next()).getSequence();
		shortest_length = tempstring.length();
		while(iter.hasNext()) {
			tempstring = ((Sequence)iter.next()).getSequence();
			if(tempstring.length() < shortest_length)
				shortest_length = tempstring.length();
		}
		return shortest_length;
	}
	
	/*
		Adds a sequence to the list of sequences
	*/
	public void addSequence(String s, String d)
	{
		Sequence seq = new Sequence(s,d);
		sequences.addLast(seq);
		size++;
		if(s.length() < length || length == 0){
			length = s.length();
		}
	}
	
	/*
	 	Finds the informative sites within this alignmentand returns them as a LinkedList<Site>
	 */
	public LinkedList<Site> findInformative()
	{
		ListIterator<Sequence> iter;
		String curr_string;
		int groupNum = code.getGroupingNumber(); // Whether the sites are comprised of single, di or tri nucleotides
		int[] recode_selection = recode_selector.getRecodeSelection(); // Array of 3 ints representing how each position in the character should be recoded
		
		for(int i = code.getOffset(); i < length; i += groupNum)
		{
			// Step 1: Create each site
			Site newSite = new Site(i+1,size);
			iter = sequences.listIterator();
			for(int j = 0; j < size; j++)
			{
				curr_string = ((Sequence)iter.next()).getSequence();
				if(i+groupNum <= length){
					String newBase = curr_string.substring(i, i+groupNum);
					newSite.addBase(newBase);
				}
			}
			
			// Step 2: Recode each site is need be
			if(groupNum == 1 && code.getSequenceType() != Code.AMINO_ACID && code.getSequenceType() != Code.STANDARD){
				newSite = code.recode(newSite,recode_selection[0],groupNum,0);
			}else if(groupNum == 2){
				// Check that we are only recoding the duplet characters
				if(newSite.getPosition() <= (length - length%2)){
					newSite = code.recode(newSite, recode_selection[0],groupNum,0);
					newSite = code.recode(newSite, recode_selection[1],groupNum,1);
				}
			}else if(groupNum == 3){
				// Check that we are only recoding the triplet characters
				if(newSite.getPosition() <= (length - length%3)){
					newSite = code.recode(newSite, recode_selection[0],groupNum,0);
					newSite = code.recode(newSite, recode_selection[1],groupNum,1);	
					newSite = code.recode(newSite, recode_selection[2],groupNum,2);	
				}
			}
			
			// Step 3: Add site to the list if it is informative and contains no unknown or missing characters. Note this still needs to be fixed up - specific for each recode and site !newSite.hasUnknown(code) && !newSite.hasMissing(code) && 
			if(newSite.getSize() > 0)
			{
				//if(!newSite.hasAmbiguous() || (recode_selection[0]+recode_selection[1]+recode_selection[2])>0){
					if(newSite.informative(code) == true)
					{	
						newSite.setChanges(code);
						informativeSites.add(newSite);
					}
					if(newSite.getChanges() == 1){
						binarySites.add(newSite);
					}
					allSites.add(newSite);
				//}
			}
		}
		return informativeSites;
	}
	
	/*
		Determines if the two sites are compatible (,1972) and returns true or false appropriately
	*/
	public boolean compatible(Site s1, Site s2)
	{
		boolean comp;
		
		// Step 1: Create list of all the pairs
		ArrayList<Pair> pairs = new ArrayList<Pair>();
		
		for(int i = 0; i < size; i++)
		{
			Pair newPair = new Pair(s1.getBase(i),s2.getBase(i));
			pairs.add(newPair);
		}
		
		// Step 2: Remove redundant pairs
		for(int i = 0; i < pairs.size(); i++)
		{
			Pair p1 = (Pair)pairs.get(i);
			int j = i+1;
			if(!p1.containsUnkOrMis(code)){
				while(j < pairs.size())
				{
					Pair p2 = (Pair)pairs.get(j);
					if(p1.equalWhole(p2) || p2.containsUnkOrMis(code))
					{
						pairs.remove(j);
					}
					else j++;
				}
			}
			else{
				pairs.remove(i);
				i--;
			}
		}
		
		// Step 3: Determine compatibility by removing pairs if either of their parts are unique to their corresponding site
		boolean unique = true;		// is changed to false if a base/amino acid matches any other in the same site
		boolean checked = false;	// used to determine whether all the most current set of pairs have been checked. Ie whether we've gone thru a cycle without removing any pairs.
		int side = 1;				// which side of the pair that is currently being analysed
		
		while(!checked)
		{
			checked = true;
			if(pairs.size()<4){
				pairs = new ArrayList<Pair>();
			}
			else{
				for(int i = 0; i < pairs.size(); i++){
					Pair p1 = (Pair)pairs.get(i);
					for(int j = 0; j < pairs.size(); j++){
						if(i!=j){
							Pair p2 = (Pair)pairs.get(j);
							if(side == 1){
								if(p1.getFirst().equals(p2.getFirst())){
									unique = false;
								}	
							}else if(side == 2){
								if(p1.getSecond().equals(p2.getSecond())){
									unique = false;
								}	
							}
						}
					}
					if(unique){
						pairs.remove(i);
						// Reset all variables to start check from the beginning
						i = 0; 
						side = 1;
						checked = false;
					}else{
						// Reset unique variable for next pair analysis
						unique = true;
						// Make sure both sides are checked before exiting while loop
						if(side==1){
							checked=false;
						}
					}
				}
				side = 3 - side; // revert to other side
			}
		}
	
		// If there are no more pairs left, site is informative
		if(pairs.size() == 0){
			comp = true;
		}else comp = false;

		return comp;
	}
	
	/*
		Informative sites are ordered based on the minimum no of changes that would have occured during 
		evolution from a common ancerstor. Changes = no of character states - 1.
		Sites are ordered from zero changes to n-1 changes, where n no of character states in the Code.
		Sites with the same no of changes are in the order they appear in the alignment.
	*/
	public LinkedList<Site> orderSitesByChanges()
	{
		ListIterator<Site> iter;
		LinkedList<Site> orderedSites = new LinkedList<Site>();
		for(int i = 1; i< size-2;i++)
		{
			iter = selectedSites.listIterator();
			while(iter.hasNext())
			{
				Site cur = iter.next();
				if(cur.getChanges() == i){
					orderedSites.add(cur);
				}
			}
		}
		return orderedSites;
	}
	
	/*
		Informative sites are ordered based on their entropy, where entropy is the freedom a site has to change.
		Sites are ordered from zero entropy (no freedom eg constant) up.
		Sites with the same entropy are in the order they appear in the alignment.
	*/
	public LinkedList<Site> orderSitesByEntropy()
	{
		ListIterator<Site> iter = selectedSites.listIterator();
		LinkedList<Site> orderedSites = new LinkedList<Site>();
		double[] ents_o = new double[informativeTotal()];
		double[] ents = new double[informativeTotal()];
		int index = 0;
		while(iter.hasNext())
		{
			Site cur = iter.next();
			ents[index] = cur.entropy(code);
			ents_o[index] = ents[index];
			index++;
		}
		Arrays.sort(ents_o);
		for(int i=0;i<ents_o.length;i++){
			iter = informativeSites.listIterator();
			index=0; 
			boolean found = false;
			while(iter.hasNext() && found==false){
				Site s = iter.next();
				if(ents_o[i]==ents[index]){
					orderedSites.add(s);
					ents[index]=-1;
					found=true;
				}
				index++;
			}
		}
		return orderedSites;
	}
	
	/* 
		Informative sites are ordered based on their overall compatability to all sites in the alignment.
		Sites are ordered from high compatibility to low compatibility.
		Sites with the same compatibility are in the order they appear in the alignment.
	*/
	public LinkedList<Site> orderSitesByComp()
	{
		calcComp();
	
		// Order sites
		ListIterator<Site> iter;
		LinkedList<Site> orderedSites = new LinkedList<Site>();
		for(int i = selectedSites.size() - 1; i >= 0 ;i--){
			iter = selectedSites.listIterator();
			while(iter.hasNext()){
				Site cur = iter.next();				
				if(cur.getHom() > i && cur.getHom() <= i+1){
					orderedSites.add(cur);
				}
			}
		}		
		return orderedSites;
	}
	
	public Map<Integer,Integer> entropyFreq()
	{	
		ListIterator<Site> iter = selectedSites.listIterator(); // maybe should be looking at all sites???
		Map<Integer,Integer> m = new HashMap<Integer,Integer>();
		while(iter.hasNext())
		{
			Site cur = iter.next();
			double ent = cur.entropy(code);
			int ent_r = (int)Math.floor(ent);
			int count = 0;
			if(m.containsKey(ent_r)){
				count = (Integer)m.get(ent_r);
			}
			m.put(ent_r,new Integer(count+1));
		}
		return m;
	}

	
	/*
		Create a plot of compatibility. 
		If passed into a BinaryGridImage will become a black and white plot. 
	*/
	public boolean[][] plot(LinkedList<Site> x, LinkedList<Site> y)
	{
		boolean[][] plot = new boolean[y.size()][x.size()];
		ListIterator<Site> x_iter;
		ListIterator<Site> y_iter = y.listIterator();
		Site y_site,x_site;
		
		for(int i = 0; y_iter.hasNext(); i++) 
		{
			//System.out.print(i + " ");
			y_site = y_iter.next();
			x_iter = x.listIterator();
			for (int j = 0; x_iter.hasNext(); j++) 
			{
				x_site = x_iter.next();
				if(j>i){
					if(compatible(y_site, x_site)) 
					{
						y_site.increaseHom();
						x_site.increaseHom();
						site_comp[i]++;
						site_comp[j]++;
						ov_comp++;
						plot[i][j] = true;
						plot[j][i] = true;
					}
					else {
						plot[i][j] = false;
						plot[j][i] = false;
					}
				}
				else if(i == j){
					plot[i][j] = true;
				}
			}
		}
		ov_comp=ov_comp/((y.size()-1)*x.size());
		COMPcalculated = true;
		return plot;
	}
	
	public void calcComp()
	{
		ListIterator<Site> iter = selectedSites.listIterator(); 
		while(iter.hasNext()){
			Site cur = iter.next();
			ListIterator<Site> iter2 = selectedSites.listIterator();
			while(iter2.hasNext()){
				Site cur2 = iter2.next();
				if(cur2.getPosition() > cur.getPosition()){
					if(compatible(cur,cur2)){
						cur.increaseHom();
						cur2.increaseHom();
					}
				}
			}
		}
		COMPcalculated = true;
	}
		
	
	public double[] NeighbourSimScore(boolean[][] plot, int gens)
	{
		double[] total_obs = new double[2];
		double[] rand = new double[gens];
		
		/*
			Calculate the observed neighbour similarity score
		*/
		double obs = 0;
		int pairs=(plot[0].length-1)*(plot[0].length-1);
		for(int i=0; i<(plot.length-1); i++){
			for(int j=i; j<plot[0].length; j++){
				if(j<(plot[0].length-1)){if(plot[i][j]==plot[i][j+1]){obs++;}}
				if(j>i){if(plot[i][j]==plot[i+1][j]){obs++;}}
			}
		}
		total_obs[0] = obs/(double)pairs;
				
		/*
			For gen Randomizations, randomize the site order and recalculate the obs similarity score.
		*/
		double max_random=0;
		int great=0;
		for(int k=0;k<gens;k++){
			double cur_obs=0;
			int[] newOrder = sample_sites(plot[0].length);
			for(int i=0; i<(plot.length-1); i++){
				for(int j=i; j<(plot[0].length); j++){
					if(j<(plot[0].length-1)){if(plot[newOrder[i]][newOrder[j]]==plot[newOrder[i]][newOrder[j+1]]){cur_obs++;}}
					if(j>i){if(plot[newOrder[i]][newOrder[j]]==plot[newOrder[i+1]][newOrder[j]]){cur_obs++;}}
				}
			}
			rand[k]=cur_obs/(double)pairs;
			if(rand[k]>max_random){max_random=rand[k];}
			if(rand[k]>total_obs[0]){great++;}
		}
		total_obs[1] = (double)great/(double)gens;
		return total_obs;		
	}

	public int[] sample_sites(int num){
		Random random = new Random();
		LinkedList<Integer> site_pos = new LinkedList();
		for(int i=0;i<num;i++){site_pos.add(new Integer(i));}
		int[] order = new int[num];
		for(int j=0;j<num;j++){
			int index = random.nextInt(site_pos.size());
			order[j] = site_pos.get(index);
			site_pos.remove(index);
		}
		return order;
	}
	
	
	/*
		Prints all the sequences to a file minus the sites that were not compatible with a 
		threshold % of other sites and have MNIC less than maxchanges.
	*/
	public void outputFile(String format, double mincompthresh, LinkedList<Double> POCthresholds, int gens, int maxchanges, String filename)
	{
		// Step 1: Determine the position of the sites in the alignment that need to be removed.
		ListIterator<Site> iter1;
		LinkedList<Integer> hetsites;
		if(POCthresholds.size()>0){
			calcPOC(gens);
		}
		else{POCthresholds.add(1.0);}
		FileOutputStream out;
		PrintStream p;
		try{
			out = new FileOutputStream(filename);
			p = new PrintStream(out);
			for(int i=0;i<POCthresholds.size();i++){
				hetsites = new LinkedList<Integer>();
				iter1 = selectedSites.listIterator();
				double POCthresh=POCthresholds.get(i);
				int cur_site = 0;
				while(iter1.hasNext())
				{
					Site s = iter1.next();
					// What proportion of other sites is this site compatible with?
					double hom = s.getHom();
					double poc = s.getPOC();
					double comp = hom/(selectedSites.size()-1);
					int changes = s.getChanges();
					
					//System.out.print(s.getPosition() + " " + s.getPOC());
					//Determining which sites (position) should be removed from the alignment based on compatibility, POC and MNIC
					if(poc > POCthresh || changes > maxchanges || comp < mincompthresh)
					{
						Integer pos = new Integer(s.getPosition());
						hetsites.addLast(pos);
						//System.out.println(" - Remove");
					}
					else{
						Integer pos = new Integer(s.getPosition());
						//System.out.println(" - Keep");
					}
					cur_site++;
				}
				
				// Step 2: Print to file all sites that do not have a position number determined to be removed.
				ListIterator<Sequence> iter2 = sequences.listIterator();
				if(format == "FASTA")
				{
					while(iter2.hasNext())
					{
						Sequence old_seq = (Sequence)iter2.next();
						String old_str = old_seq.getSequence();
						String new_str = "";
						int start = 0;
						ListIterator iter3 = hetsites.listIterator();
						while(iter3.hasNext())
						{
							Integer site = (Integer)iter3.next();
							new_str+=old_str.substring(start,site-1);
							start = site;
						}
						new_str+=old_str.substring(start);
						p.print(">");
						p.println(old_seq.getDescription());
						p.println(new_str);
					}
					p.println("");
				}
				else if(format == "SPHYLIP")
				{
					p.println(String.format(" %d %d",size,(length-hetsites.size())));
					while(iter2.hasNext())
					{
						Sequence old_seq = (Sequence)iter2.next();
						String old_str = old_seq.getSequence();
						String new_str = "";
						int start = 0;
						ListIterator iter3 = hetsites.listIterator();
						while(iter3.hasNext())
						{
							Integer site = (Integer)iter3.next();
							new_str+=old_str.substring(start,site-1);
							start = site;
						}
						new_str+=old_str.substring(start);
						String descrip = old_seq.getDescription();
						String blanks = "          ";
						if(descrip.length()>9){
							descrip = String.format("%s%s",descrip.substring(0,9)," ");
						}
						else if(descrip.length()<10){
							int n = 10-descrip.length();
							blanks=blanks.substring(0,n);
							descrip=String.format("%s%s",descrip,blanks);
						}
						p.print(descrip);
						p.println(new_str);
					}
				}
				else if(format == "NEXUS"){}
				else if(format == "IPHYLIP"){}
			}p.close();
		}catch(IOException ioe){System.err.println("Error writing to file");}
	}
	
	/*
		Print out all the informative sites in the alignment
	*/
	public void printInformative()
	{
		for(int i = 0; i < informativeSites.size(); i++)
		{
			Site aSite = (Site)informativeSites.get(i);
			System.out.print(aSite.getPosition());
		}
	}
		
	public void calcPOC(int reps)
	{
		LinkedList<double[]> comp_list =  new LinkedList();
		double[] POC = new double[selectedSites.size()];
		ListIterator<Site> iter = selectedSites.listIterator();
		while(iter.hasNext()){
			Site cur = iter.next();
			double[] comp_prop = new double[reps];
			int count = 0;
			double cur_poc = 0;
			for(int k=1;k<(reps+1);k++){
				String[] neworder = cur.randomize();
				Site rand = new Site(cur.getPosition(),size);
				rand.addMembers(neworder);
				ListIterator<Site> iter2 = selectedSites.listIterator();
				int rand_comp=0;
				while(iter2.hasNext()){
					Site cur_site=iter2.next();
					if(rand.getPosition() != cur_site.getPosition()){
						if(compatible(rand,cur_site)){
							rand_comp++;
						}
					}
				}
				if(cur.getHom() <= rand_comp){
					cur_poc++;
				}
			}
			double temp = cur_poc/(double)reps;
			cur.setPOC(temp);
		}
		POCcalculated = true;
	}
	
	/*
		Calculates the overall compatibility of each site within the alignment
	*/
	public double ovSiteComp(Site s)
	{
		double count = 0;
		ListIterator iter = uniqueSites.listIterator();
		while(iter.hasNext())
		{
			Site cur = (Site)iter.next();
			if(s.getPosition() != cur.getPosition()){
				if(compatible(s,cur)){
					count+=uniquemap.get(cur.getUniquePos());
				}
			}
		}
		count+=uniquemap.get(s.getUniquePos())-1;
		count=count/(selectedSites.size()-1);
		return count;
	}
	
	/*
		Helper method the round down each sites significance for printing
	*/
	public double roundSig(double sig,int dp)
	{
		double mult=Math.pow(10,dp);
		return ((double)((int)Math.round(mult*sig)))/(int)mult;
	}
	
	/*
		Determines whether two sites are the same in regards to character composition and position
	*/
	public boolean sameSite(Site s1, Site s2)
	{
		int mnic1 = s1.getChanges();
		int mnic2 = s2.getChanges();
		boolean same=true;
		if(mnic1!=mnic2){same = false;}
		else{
			int[] comps1= s1.getSplits(code);
			int[] comps2= s2.getSplits(code);
			if(comps1.length==comps2.length){
				for(int i=0;i<comps1.length;i++){
					if(comps1[i]!=comps2[i]){same = false;}
				}
				if(same){
					LinkedList<Pair> pairs = new LinkedList<Pair>();
					for(int i = 0; i < size; i++){
						Pair newPair = new Pair(s1.getBase(i),s2.getBase(i));
						pairs.add(newPair);
					}
					ListIterator iter = pairs.listIterator();
					for(int i = 0; i < pairs.size(); i++){
						Pair p1 = (Pair)pairs.get(i);
						int j = i+1;
						if(!p1.containsUnkOrMis(code)){
							while(j < pairs.size()){
								Pair p2 = (Pair)pairs.get(j);
								if(p1.equalWhole(p2) || p2.containsUnkOrMis(code)){pairs.remove(j);}
								else j++;
							}
						}
						else{
							pairs.remove(i);
							i--;
						}
					}
					if(pairs.size()!=(mnic1+1)){
						same=false;
					}
				}
			}else{same=false;}
		}
		return same;
	}
	
	public void setSelected(LinkedList<Site> list){
		selectedSites = list;
	}
	
	/*
		Finds all the unique informative sites ie removes redundant sites with the same character composition and pattern
	*/
	public void findUnique()
	{
		ListIterator iter = selectedSites.listIterator();
		uniqueSites.add((Site)iter.next());
		uniquemap.add(new Integer(1));
		while(iter.hasNext()){
			Site s1 = (Site)iter.next();
			ListIterator iter2 = uniqueSites.listIterator();
			boolean found = false;
			int pos=0;
			while(iter2.hasNext() && !found){
				Site s2 = (Site)iter2.next();
				if(sameSite(s1,s2)){
					found=true;
					uniquemap.set(pos,(uniquemap.get(pos)+1));
					s2.setUniquePos(pos);
				}
				pos++;
			}
			if(!found){
				uniqueSites.add(s1);
				uniquemap.add(new Integer(1));}
		}
	}
	
	public void printSiteSummary()
	{
		FileOutputStream out;
		PrintStream p;
		try{
			out = new FileOutputStream("site_summary.csv");
			p = new PrintStream(out);
			//if(!POCcalculated){calcPOC(1000);}
			if(!COMPcalculated){calcComp();}
			ListIterator<Site> iter = selectedSites.listIterator();
			p.println("Site No.,MNIC,Compatibility,POC,Entropy"); 
			while(iter.hasNext()){
				Site cur = iter.next();
				p.println(cur.getPosition() + "," + cur.getChanges() + "," + (double)cur.getHom()/(selectedSites.size()-1) + "," + cur.getPOC() + "," +cur.entropy(code));
			}
		}
		catch(IOException ioe)
		{
			System.err.println("Error writing to file");
		}
	}
		
}
