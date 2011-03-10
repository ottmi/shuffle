/*
	Class which defines a sequence - a string of nucleotide or amino acids 
*/

import java.util.*;
import java.io.*;

public class Sequence
{
	private String seq; // string of nucleotide or amino acids 
	private String description; // Eg. species and gene to which the sequence belongs
	
	public Sequence(String s, String d)
	{
		seq = s;
		description = d;
	}
	
	/*
		Simple 'get' methods
	*/
	public String getDescription(){return description;}
	public String getSequence(){return seq;}
	
	/*
		Resets the seq string
	*/
	public void setSequence(String sequence) { seq = sequence; }
}
		