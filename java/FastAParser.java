/*
	A class to read in sequence data adhering to the FastA format, as defined
	at http://www.ncbi.nlm.nih.gov/blast/fasta.shtml
*/

import java.io.*;
import java.util.*;


public class FastAParser {
	
	LinkedList<Sequence> sequences;
	BufferedReader input;
	
	public FastAParser(File file) {
		sequences = new LinkedList<Sequence>();
		try {
			input = new BufferedReader(new FileReader(file));
		}
		catch(FileNotFoundException fnfe) {
			System.err.println("file "+file.getName()+" not located, exiting");
			fnfe.printStackTrace();
			System.exit(1);
		}
	}
	public FastAParser(String filename) {
		sequences = new LinkedList<Sequence>();
		
		//instantiates the BufferedReader
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
		
		//moves onto first line of file
		try {
			currstr = input.readLine();
		}
		catch(IOException ioe) {
			System.err.println("ERROR encountered reading in line, exiting");
			ioe.printStackTrace();
			System.exit(1);
		}
		
		//Main read-in loop
		String tempstr = new String();
		boolean finished = false;
		String descrip = currstr.substring(1);
		while(!finished) {
			//moves onto next line of input file
			try {
				currstr = input.readLine();
			}
			catch(IOException ioe) {
				System.err.println("ERROR encountered reading in line, exiting");
				ioe.printStackTrace();
				System.exit(1);
			}
			
			//checks to see if string is empty, if so, end of file, add tempstr to sequences, break
			if(currstr == null) {
				Sequence seq = new Sequence(tempstr,descrip);
				sequences.addLast(seq);
				finished = true;
				continue;
			}
			
			//if first character is '>', add tempstr to sequences list, wipe tempstr, then break loop
			if(currstr.charAt(0) == '>') {
				Sequence seq = new Sequence(tempstr,descrip);
				sequences.addLast(seq);
				descrip = currstr.substring(1);
				tempstr = new String();
				continue;
			}
			
			//attach currstr input string to end of tempstr
			tempstr += removeSpace(currstr.toUpperCase()); //FA convention
		}
		return sequences;
	}
		
	/*
		Method to remove any white space contained within the sequence.
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