/*
	A class to read in sequence data adhering to the Nexus format, as defined
	at http://www.ncbi.nlm.nih.gov/blast/fasta.shtml
*/

import java.util.*;
import java.io.*;
import java.lang.*;

public class NexusParser
{
	LinkedList<Sequence> sequences;
	BufferedReader input;
	
	public NexusParser(File file)
	{
		sequences = new LinkedList<Sequence>();
		
		try {
			input = new BufferedReader(new FileReader(file));
		}
		catch(FileNotFoundException fnfe) {
			System.err.println("file " + file.getName() + " not located, exiting");
			fnfe.printStackTrace();
			System.exit(1);
		}
	}
	
	public NexusParser(String filename) 
	{
		sequences = new LinkedList<Sequence>();
		
		try {
			input = new BufferedReader(new FileReader(new File(filename)));
		}
		catch(FileNotFoundException fnfe) {
			System.err.println("file "+filename+" not located, exiting");
			fnfe.printStackTrace();
			System.exit(1);
		}
	}
	
	/*
		Parses the file, collecting and storing the sequence and description in a new Sequence.
	*/
	public LinkedList<Sequence> parseFile() 
	{
		String currstr = null;
		
		// moves onto first line of file
		try {
			currstr = input.readLine();
		}
		catch(IOException ioe) {
			System.err.println("ERROR encountered reading in line, exiting");
			ioe.printStackTrace();
			System.exit(1);
		}
		if(currstr.length() > 0){
				currstr = currstr.toUpperCase();
		}
		
		boolean matrixFound = false;
		int taxa = 0;
		char unknown = '?';
		char gap = '-';
		boolean interleaved = false;
		
		// Step 1: collect all necessary data from the file, including the number of taxa,
		// the symbols for a missing characters and a gap and whether the file is interleaved or not
		// This is done until the word 'Matrix' is read in - Matrix indicates the start of the sequences.
		while(!matrixFound)
		{
			if(currstr.contains("NTAX"))
			{
				int ind1a = currstr.indexOf("NTAX");
				String temp = currstr.substring(ind1a+5);
				int ind1b = temp.indexOf(";");
				int ind1c = temp.indexOf(" ");
				if(ind1c < ind1b && ind1c > 0)
					taxa = new Integer(temp.substring(0,ind1c));
				else taxa = new Integer(temp.substring(0,ind1b));
			}
			if(currstr.contains("MATRIX"))
			{
				matrixFound = true;
			}
			if(currstr.contains("MISSING="))
			{
				int ind2 = currstr.indexOf("MISSING");
				unknown = currstr.charAt(ind2+8);
			}
			if(currstr.contains("GAP="))
			{
				int ind3 = currstr.indexOf("GAP");
				gap = currstr.charAt(ind3+4);
			}
			if(currstr.contains("INTERLEAVE"))
			{
				int ind4 = currstr.indexOf("INTERLEAVE");
				if(currstr.charAt(ind4+11) != 'N'){
					interleaved = true;
				}
			}
			
			// read in next line of file
			try {
				currstr = input.readLine();
			}
			catch(IOException ioe) {
				System.err.println("ERROR encountered reading in line, exiting");
				ioe.printStackTrace();
				System.exit(1);
			}
			if(currstr.length() > 0){
				currstr = currstr.toUpperCase();
			}
		}
		
		// Step 2: Matrix is found, so begin reading in sequences.
		int counttaxa = 0;
		boolean finished = false;
		String tempstr = new String();
		String descrip = new String();
		ArrayList<String> lines = new ArrayList<String>();
	
		// All non-blank lines are read in and stored in a list - lines.
		// Lines can not contain the character '[' as this is sometimes used to indicate how many characters are in that particular block
		while(!finished)
		{

			if(currstr.length() > 1 && currstr.charAt(0) != '['){
				currstr = currstr.toUpperCase();
				if(currstr.contains(";") || currstr == null){
					finished = true;
				}
				else{
					lines.add(currstr);
				}
			}
			if(!finished)
			{
				try {
					currstr = input.readLine();
				}
				catch(IOException ioe) {
					System.err.println("ERROR encountered reading in line, exiting");
					ioe.printStackTrace();
					System.exit(1);
				}
			}
		}
		
		// Step 3: Determine how many lines each sequence covers
		int numlines = lines.size()/taxa;
		
		// If the file is interleaved
		if(interleaved)
		{
			// Step 4a: For each taxa, concat each of the lines together to form the sequence.
			// For example, if taxa1 had the first sequence in the file and there where 10 taxa in the file, concat lines 1,11,21...numlines*10+1
			while(counttaxa < taxa)
			{
				for(int i = 0; i< numlines; i++)
				{
					String line = lines.get(taxa*i + counttaxa);
					String tempshort;
					int start = line.indexOf(" ");
					if(line.indexOf("	") < start && line.indexOf("	") > 0){
						start = line.indexOf("	");
					}
					if(i==0)
					{
						descrip = line.substring(0,start);
					}
					if(line.charAt(line.length()-1) == '\n'){
						tempshort = removeSpace(line.substring(start+1, line.length()-1));
					}
					else {tempshort = removeSpace(line.substring(start+1, line.length()));}
					tempstr += tempshort;
				}
				Sequence s = new Sequence(tempstr,descrip);
				sequences.add(s);
				tempstr = "";
				counttaxa++;
			}
		}
		// If the file is sequential
		else
		{
			// Step 4b: For each taxa, concat each of the lines together to form the sequence.
			while(counttaxa < taxa)
			{	
				
				for(int i = 0; i< numlines; i++)
				{	
					String line = lines.get(counttaxa*(numlines) + i);
					int start = line.indexOf(" ");
					if(i == 0)
					{
						descrip = line;
					}
					else
					{
						String tempshort = removeSpace(line);
						tempstr += tempshort;
					}
				}
				Sequence s = new Sequence(tempstr,descrip);
				sequences.add(s);
				tempstr = "";
				counttaxa++;
			}
		}
		return sequences;
	}
	
	/*
		Method to remove any white space contained within the sequence.
		Eg. If characters in sequence are grouped in tens, with a space separating the groups, this method will remove the spaces.
	*/
	public String removeSpace(String oldString)
	{
		String str = new String();
		String[] newString = oldString.split(" ");
		for(int i = 0; i <  newString.length; i++)
		{
			str += newString[i];
		}
		return str;
	}
}