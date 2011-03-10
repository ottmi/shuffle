/*
	Pair is a helper class used for the analysis of compatibility between sites.
	A Pair is made up of two characters from a single sequence corresponding to two sites being analysed for compatibility
	In an alignment of n sequences, with m sites for any compatibility analysis between two sites
	there are n Pairs.
*/

import java.util.*;

public class Pair 
{
	private String first; // 1st member of the pair
	private String second; // 2nd member of the pair
	
	public Pair(String a, String b)
	{
		first = a;
		second = b;
	}
	
	public String getFirst()
	{
		return first;
	}
	
	public String getSecond()
	{
		return second;
	}
	
	/*
		Determines if the Pair is equal to another pair.
		Used in the Alignment:compatible() method to remove redundant pairs
	*/
	public boolean equalWhole(Pair p2)
	{
		if(first.equals(p2.first) && second.equals(p2.second))
		{
			return true;
		}
		else return false;
	}
	
	/*
		Determines if the Pair contains an unknown or missing character.
		Used in the Alignment:findInformative() method
	*/
	public boolean containsUnkOrMis(Code code)
	{
		boolean conts = false;
		if(first.contains(String.valueOf(code.getUnknown())) || first.contains(String.valueOf(code.getMissing()))){
			conts = true;
		}
		else if(second.contains(String.valueOf(code.getUnknown())) || second.contains(String.valueOf(code.getMissing()))){
			conts = true;
		}
		else if(code.getSequenceType()>0){
			if(containsAmbig()){
				conts = true;
			}
		}
		return conts;
	}
	
	/*
		For nucleotide alignments, determines if the Pair contains an ambiguous character.
		Used in the Alignment:findInformative() method as Shuffle regards ambigous characters as unknown
	*/
	public boolean containsAmbig()
	{
		char[] ambig_chars = {'R','Y','K','M','S','W','B','D','H','V'};
		boolean ambig=false;
		for(int j =0;j<ambig_chars.length;j++){
			if(first.contains(String.valueOf(ambig_chars[j])) || second.contains(String.valueOf(ambig_chars[j])))
			{
				ambig = true;
			}
		}
		return ambig;
	}
}
