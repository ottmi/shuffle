/*
	Class which represents a site in an alignment.
	By definition a site is a row in an alignment - ie. all the characters at a single position in an alignment
*/

import java.util.*;
import java.io.*;
import java.lang.*;

public class Site 
{
	private int number; // ID number for the site (position within the alignment)
	public int size; // number of characters making up the site
	private String[] members; // characters making up the site
	private int hom; // no of sites this site is compatible with in the alignment
	private int changes; // the MNIC expressed by this site
	private int completeness; // number of characters within the site that are known (not unknown, missing or ambiguous)
	private double POC; // the significance of the site's compatibility score
	private int unique_pos;
	
	public Site(int num, int seqs)
	{
		number = num;
		size = 0;
		members = new String[seqs];
		changes = 0;
		hom = 0;
		completeness=0;
		POC=0;
		unique_pos=0;
	}
	
	/*
		Adds a new character to the Site
		Called in the Alignment:findInformative() method
	*/
	public void addBase(String c)
	{
		members[size] = c;
		size++;
	}
	
	/*
		Some get methods
	*/
	public String getBase(int i){return members[i];}
	public int getPosition(){return number;}
	public int getHom(){return hom;}
	public int getChanges(){return changes;}
	public int getSize(){return size;}
	public double getPOC(){return POC;}
	public int getCompleteness(){ return completeness;}
	public double getOvComp(int sitetot){return ((double)hom/(sitetot-1))*100;}
	public void setUniquePos(int pos){unique_pos = pos;}
	public int getUniquePos(){return unique_pos;}
	
	/*
		Internal counter for Alignment methods, orderSitesByCompatibility() and outputFile() 
	*/
	public void increaseHom(){hom++;}
	
	/*
		Checks if the site is parsimony informative.
		A site is defined as parsimony informative is it has two or more character groups with two or more members
	*/
	public boolean informative(Code code)
	{
		setChanges(code);
		calcCompleteness2(code);
		int count = 0;
		int i = 0;
		String b = null;
		while(count < 2 && i < size)
		{
			int j = i + 1; 
			boolean equal = false;
			while(j < size && equal == false)
			{
				if(members[i].equals(members[j])){
					if(!members[j].contains(String.valueOf(code.getUnknown())) && !members[j].contains(String.valueOf(code.getMissing()))
							&& !members[i].contains(String.valueOf(code.getUnknown())) && !members[i].contains(String.valueOf(code.getMissing())))
					{
						if(code.getSequenceType() == 0 || code.getSequenceType() > 4 || (!isAmbiguous(members[i]) && !isAmbiguous(members[j]))){
							if(b == null){
								b = members[i];
								count++;
							}
							else if(!(b.equals(members[i]))){
								count++;
							}
							equal = true;
						}
						else j++;
					}else j++;
				}else j++;
			}
			i++;
		}
		if(count == 2){
			return true;
		}
		else return false;
	}
	
	/*
		Sets the value for the variable changes by comparing each member of the site to each member of the Code.
		changes ~ MNIC i.e., the number of different characters minus 1
	*/
	public void setChanges(Code aCode)
	{
		String[] mems = aCode.getMembers2(); //creates a linkedlist containing all the possible characters that might appear
		changes = -1; 
		
		for(int mems_index = 0; mems_index < mems.length; mems_index++)
		{
			//iterates through the array of member strings of this site, if a match is found, changes++
			for(int i = 0; i < members.length; i++) 
			{
				if((members[i].equals(mems[mems_index]))) 
				{
					changes++;
					i = members.length;
				}
			}
		}
	}
	
	/*
		Determines whether the site contains an unknown character.
	*/
	public boolean hasUnknown(Code aCode)
	{
		boolean unknown = false;
		char unk = aCode.getUnknown();
		String unknown_char = String.valueOf(unk);
		
		for(int i = 0; i<size; i++)
		{
			if(members[i].contains(unknown_char))
			{
				unknown = true;
			}
		}
		return unknown;
	}
	
	/*
		Determines whether the site contains a missing character.
	*/
	public boolean hasMissing(Code aCode)
	{
		boolean missing = false;
		char mis = aCode.getMissing();
		String missing_char = String.valueOf(mis);
		
		for(int i = 0; i<size; i++)
		{
			if(members[i].contains(missing_char))
			{
				missing = true;
			}
		}
		return missing;
	}
	
	/*
		Determines whether the site contains an ambiguous character. Only for DNA.
	*/
	public boolean hasAmbiguous()
	{
		char[] ambig_chars = {'R','Y','K','M','S','W','B','D','H','V'};
		boolean ambig=false;
		for(int i = 0; i<size; i++)
		{
			for(int j =0;j<ambig_chars.length;j++){
				if(members[i].contains(String.valueOf(ambig_chars[j])))
				{
					ambig = true;
				}
			}
		}
		return ambig;
	}
	
	/*
		Determines whether a particular character is ambiguous. Only for DNA.
	*/
	public boolean isAmbiguous(String s)
	{
		char[] ambig_chars = {'R','Y','K','M','S','W','B','D','H','V'};
		boolean ambig=false;
		for(int i =0;i<ambig_chars.length;i++){
			if(s.contains(String.valueOf(ambig_chars[i])))
			{
				ambig = true;
				i=ambig_chars.length;
			}
		}
		return ambig;
	}
	
	/*
		Converts the array of characters making up the site to a string.
	*/
	public String toString() 
	{
		String rString = new String();
		for(int i = 0; i < members.length; i++) {
			rString += members[i];
		}
		return rString;
	}
	
	/*
		Calculates the entropy of a site
		The entropy of a site is a measure of its freedom to change.
		Formula used = ln(n!/(prod(r!) x prod(t!))), 
		where:
		n is the number of sequences making up the site, 
		r is a vector with the frequency of unique characters at the site,
		t is a vector with the frequency of unique frequencies in r
		For example, Site x = ACGTGTAT, n = 8, r = [2,1,2,3], t = [1,2,1]
	*/
	public double entropy(Code code)
	{
		System.out.print("Calculating Entropy of ");
		print();
		double ent = 0;
		int s = size; // number of characters making up the site
		int[] r = getSplits(code); // vector with the frequency of unique characters
		int max = findMax(r);
		double denom=1;
		Map<Integer,Integer> m = new HashMap<Integer,Integer>();
		// loop goes through r vector counting the frequency of each unique value
		for(int i=0;i<r.length;i++){
			int count=0;
			if(m.containsKey(r[i])){
				count=(Integer)m.get(r[i]);
			}
			m.put(r[i],new Integer(count+1));
			denom*=factorial(r[i],0);
			
		}
		int[] t = new int[m.size()];
		int index=0;
		Iterator iter = m.entrySet().iterator();
		while(iter.hasNext())
		{
			Map.Entry pair = (Map.Entry)iter.next();
			Integer value = (Integer)pair.getValue();
			t[index]=value;
			denom*=factorial(value,0);
			index++;
		}
		for(int j=0;j<r.length;j++){

			System.out.print(r[j] + ",");
		}
		System.out.println("");
		for(int k=0;k<t.length;k++){
			System.out.print(t[k] + ",");
		}
		//ent=Math.log(factorial(completeness,0)/denom);
		ent=Math.log(factorial(s,0)/denom);
		System.out.println("");
		System.out.println(ent);
		return ent;
	}
	
	/*
		Helper method for entropy()
	*/
	public double factorial(int n, int end)
	{
		double result = 1;
		for(int i=n;i>end;i--){
			result=result*i;
		}
		return result;
	}
	
	/*
		Helper method for entropy()
	*/
	public int findMax(int[] comps)
	{	
		int max = 0;
		for(int i=0;i<comps.length;i++)
		{
			if(comps[i]>max){
				max = comps[i];
			}
		}
		return max;
	}

	
		/*
		Calculates the completeness of the site ie. what proportion of sites are not missing, unknown or ambiguous
	*/
	public void calcCompleteness2(Code code)
	{
		int count = size;
		for(int i=0;i<members.length;i++){
			if(members[i].contains(String.valueOf(code.getUnknown())) || members[i].contains(String.valueOf(code.getMissing()))){
				count--;
			}
			/*else if(code.getSequenceType()>0 && code.getSequenceType()<4){
				for(int j=0;j<(code.NUCLEOTIDE_CHARACTERS_SINGLE_AMBIGUOUS).length; j++){
					if(members[i].contains(String.valueOf(code.NUCLEOTIDE_CHARACTERS_SINGLE_AMBIGUOUS[j]))){
						count--;
						break;
					}
				}
			}*/
		}
		completeness=count;
	}
	
	/*
		Determines the character splits of the site. eg [A,A,A,G,G,A,A,T,C] would have splits (5|2|1|1)
		Helper method for entropy()
	*/
	public int[] getSplits(Code code)
	{
		Map<String,Integer> m= new HashMap<String,Integer>();
		for(int i=0;i<members.length;i++){
			int count = 0;
			//if(!members[i].contains(String.valueOf(code.getUnknown())) && !members[i].contains(String.valueOf(code.getMissing()))){
				if(m.containsKey(members[i])){
					count=(Integer)m.get(members[i]);
				}
				m.put(members[i],new Integer(count+1));
			//}
		}
		int[] comps = new int[m.size()];
		int index=0;
		Iterator iter = m.entrySet().iterator();
		while(iter.hasNext())
		{
			Map.Entry pair = (Map.Entry)iter.next();
			Integer value = (Integer)pair.getValue();
			comps[index] = value;
			index++;
		}
		Arrays.sort(comps);
		return comps;
	}

	
	public boolean isComplete()
	{
		if(completeness < size){
			return false;
		}else{ return true;}
	}
	
	/*
		Randomizes the order of known characters within the site
		Used in calculating the probability of compatibility of the site in the alignment
	*/
	public String[] randomize()
	{
		Random gen = new Random();
		LinkedList<Integer> positions = new LinkedList();
		for(int i=0;i<members.length;i++){
			positions.add(i);
		}
		String[] reorder = new String[members.length];
		for(int j=0;j<members.length;j++){
			int pos = gen.nextInt(positions.size());
			reorder[j] = members[positions.get(pos)];
			positions.remove(pos);
		}
		return reorder;
	}
	
	/*
		Adds/changes members of a site
		Used in the randomization process
	*/
	public void addMembers(String[] mems){
		members = mems;
	}
	
	/*
		Sets the probability of compatibility (significance) for this site
	*/
	public void setPOC(double poc){
		POC=poc;
	}
	
	/*
		Calculates codon position. Nucleotides triplets Only!
	*/
	
	public int getCP(int datatype)
	{
		int cp = number%datatype;
		if(cp == 0){cp=datatype;}
		return cp;
	}
	
	/*
		Prints the Site to standard out
	*/
	public void print(){
		System.out.println(toString());
	}
}



