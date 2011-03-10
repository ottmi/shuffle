/*
	Recode Selector is the GUI panel where users decide whether they want to recode sites and how.
	There are 13 different ways nucleotides can be recoded.
	If a site is made up of more than one character, each position within the site can be recoded separately 
*/

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class RecodeSelector extends JPanel implements ActionListener 
{
	// Static recoding variables
	public static final int GROUP_NONE = 0;
	public static final int GROUP_PURINES = 1;
	public static final int GROUP_PYRIMIDINES = 2;
	
	public static final int GROUP_PURINES_AND_PYRIMIDINES = 3;
	
	public static final int GROUP_KETONES = 4;
	public static final int GROUP_AMINOS = 5;
	
	public static final int GROUP_KETONES_AND_AMINOS = 6;
	
	public static final int GROUP_STRONG = 7;
	public static final int GROUP_WEAK = 8;
	
	public static final int GROUP_STRONG_AND_WEAK = 9;
	
	public static final int GROUP_NOT_A = 10;
	public static final int GROUP_NOT_C = 11;
	public static final int GROUP_NOT_G = 12;
	public static final int GROUP_NOT_T = 13;
	
	JRadioButton[] gp1; // radiobuttons for position 1
	JRadioButton[] gp2; // radiobuttons for position 2
	JRadioButton[] gp3; // radiobuttons for position 3
	
	private int[] current_selection; // what the radio buttons currently show as the recoding scheme
	
	public RecodeSelector() {
		super(new BorderLayout());
		current_selection = new int[3];
		current_selection[0] = RecodeSelector.GROUP_NONE;
		current_selection[1] = RecodeSelector.GROUP_NONE;
		current_selection[2] = RecodeSelector.GROUP_NONE;
		
		/*
			creates radiobuttons with 3 lots of each possible scheme allowing seperate recoding of the 
			max 3 positions a site can have
		*/
		
		// Group 1 = Position 1
		gp1 = new JRadioButton[14];
		
		JRadioButton group_none_button = new JRadioButton("A,C,G,T");
		group_none_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_NONE));
		group_none_button.setSelected(true);
		gp1[0] = group_none_button;
		
		JRadioButton group_purines_button = new JRadioButton("R,C,T");
		group_purines_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_PURINES));
		gp1[1] = group_purines_button;
		
		JRadioButton group_pyrimidines_button = new JRadioButton("Y,A,G");
		group_pyrimidines_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_PYRIMIDINES));
		gp1[2] = group_pyrimidines_button;
		
		JRadioButton group_purines_and_pyrimidines_button = new JRadioButton("R,Y");
		group_purines_and_pyrimidines_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_PURINES_AND_PYRIMIDINES));
		gp1[3] = group_purines_and_pyrimidines_button;
		
		JRadioButton group_ketones_button = new JRadioButton("K,A,C");
		group_ketones_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_KETONES));
		gp1[4] = group_ketones_button;
		
		JRadioButton group_aminos_button = new JRadioButton("M,G,T");
		group_aminos_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_AMINOS));
		gp1[5] = group_aminos_button;
		
		JRadioButton group_ketones_and_aminos_button = new JRadioButton("K,M");
		group_ketones_and_aminos_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_KETONES_AND_AMINOS));
		gp1[6] = group_ketones_and_aminos_button;
		
		JRadioButton group_strong_button = new JRadioButton("S,A,T");
		group_strong_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_STRONG));
		gp1[7] = group_strong_button;
		
		JRadioButton group_weak_button = new JRadioButton("W,C,G");
		group_weak_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_WEAK));
		gp1[8] = group_weak_button;
		
		JRadioButton group_strong_and_weak_button = new JRadioButton("S,W");
		group_strong_and_weak_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_STRONG_AND_WEAK));
		gp1[9] = group_strong_and_weak_button;
		
		JRadioButton group_not_a_button = new JRadioButton("B,A");
		group_not_a_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_A));
		gp1[10] = group_not_a_button;
		
		JRadioButton group_not_c_button = new JRadioButton("D,C");
		group_not_c_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_C));
		gp1[11] = group_not_c_button;
		
		JRadioButton group_not_g_button = new JRadioButton("H,G");
		group_not_g_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_G));
		gp1[12] =group_not_g_button;
		
		JRadioButton group_not_t_button = new JRadioButton("V,T");
		group_not_t_button.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_T));
		gp1[13] = group_not_t_button;
		
		// Group 2 = Position 2
		gp2 = new JRadioButton[14];
		JRadioButton group_none_button2 = new JRadioButton("A,C,G,T");
		group_none_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_NONE+14));
		group_none_button2.setSelected(true);
		gp2[0] = group_none_button2;
		
		JRadioButton group_purines_button2 = new JRadioButton("R,C,T");
		group_purines_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_PURINES+14));
		gp2[1] = group_purines_button2;
		
		JRadioButton group_pyrimidines_button2 = new JRadioButton("Y,A,G");
		group_pyrimidines_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_PYRIMIDINES+14));
		gp2[2] = group_pyrimidines_button2;
		
		JRadioButton group_purines_and_pyrimidines_button2 = new JRadioButton("R,Y");
		group_purines_and_pyrimidines_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_PURINES_AND_PYRIMIDINES+14));
		gp2[3] = group_purines_and_pyrimidines_button2;
		
		JRadioButton group_ketones_button2 = new JRadioButton("K,A,C");
		group_ketones_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_KETONES+14));
		gp2[4] = group_ketones_button2;
		
		JRadioButton group_aminos_button2 = new JRadioButton("M,G,T");
		group_aminos_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_AMINOS+14));
		gp2[5] = group_aminos_button2;
		
		JRadioButton group_ketones_and_aminos_button2 = new JRadioButton("K,M");
		group_ketones_and_aminos_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_KETONES_AND_AMINOS+14));
		gp2[6] = group_ketones_and_aminos_button2;
		
		JRadioButton group_strong_button2 = new JRadioButton("S,A,T");
		group_strong_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_STRONG+14));
		gp2[7] = group_strong_button2;
		
		JRadioButton group_weak_button2 = new JRadioButton("W,C,G");
		group_weak_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_WEAK+14));
		gp2[8] = group_weak_button2;
		
		JRadioButton group_strong_and_weak_button2 = new JRadioButton("S,W");
		group_strong_and_weak_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_STRONG_AND_WEAK+14));
		gp2[9] = group_strong_and_weak_button2;
		
		JRadioButton group_not_a_button2 = new JRadioButton("B,A");
		group_not_a_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_A+14));
		gp2[10] = group_not_a_button2;
		
		JRadioButton group_not_c_button2 = new JRadioButton("D,C");
		group_not_c_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_C+14));
		gp2[11] = group_not_c_button2;
		
		JRadioButton group_not_g_button2 = new JRadioButton("H,G");
		group_not_g_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_G+14));
		gp2[12] = group_not_g_button2;
		
		JRadioButton group_not_t_button2 = new JRadioButton("V,T");
		group_not_t_button2.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_T+14));
		gp2[13] = group_not_t_button2;
		
		// Group 3 = Position 3
		gp3 = new JRadioButton[14];
		JRadioButton group_none_button3 = new JRadioButton("A,C,G,T");
		group_none_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_NONE+28));
		group_none_button3.setSelected(true);
		gp3[0] = group_none_button3;
		
		JRadioButton group_purines_button3 = new JRadioButton("R,C,T");
		group_purines_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_PURINES+28));
		gp3[1] = group_purines_button3;
		
		JRadioButton group_pyrimidines_button3 = new JRadioButton("Y,A,G");
		group_pyrimidines_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_PYRIMIDINES+28));
		gp3[2] = group_pyrimidines_button3;
		
		JRadioButton group_purines_and_pyrimidines_button3 = new JRadioButton("R,Y");
		group_purines_and_pyrimidines_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_PURINES_AND_PYRIMIDINES+28));
		gp3[3] = group_purines_and_pyrimidines_button3;
		
		JRadioButton group_ketones_button3 = new JRadioButton("K,A,C");
		group_ketones_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_KETONES+28));
		gp3[4] = group_ketones_button3;
		
		JRadioButton group_aminos_button3 = new JRadioButton("M,G,T");
		group_aminos_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_AMINOS+28));
		gp3[5] = group_aminos_button3;
		
		JRadioButton group_ketones_and_aminos_button3 = new JRadioButton("K,M");
		group_ketones_and_aminos_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_KETONES_AND_AMINOS+28));
		gp3[6] = group_ketones_and_aminos_button3;
		
		JRadioButton group_strong_button3 = new JRadioButton("S,A,T");
		group_strong_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_STRONG+28));
		gp3[7] = group_strong_button3;
		
		JRadioButton group_weak_button3 = new JRadioButton("W,C,G");
		group_weak_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_WEAK+28));
		gp3[8] = group_weak_button3;
		
		JRadioButton group_strong_and_weak_button3 = new JRadioButton("S,W");
		group_strong_and_weak_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_STRONG_AND_WEAK+28));
		gp3[9] = group_strong_and_weak_button3;
		
		JRadioButton group_not_a_button3 = new JRadioButton("B,A");
		group_not_a_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_A+28));
		gp3[10] = group_not_a_button3;
		
		JRadioButton group_not_c_button3 = new JRadioButton("D,C");
		group_not_c_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_C+28));
		gp3[11] = group_not_c_button3;
		
		JRadioButton group_not_g_button3 = new JRadioButton("H,G");
		group_not_g_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_G+28));
		gp3[12] = group_not_g_button3;
		
		JRadioButton group_not_t_button3 = new JRadioButton("V,T");
		group_not_t_button3.setActionCommand(Integer.toString(RecodeSelector.GROUP_NOT_T+28));
		gp3[13] = group_not_t_button3;
		
		// Puts the buttons into a ButtonGroup
		ButtonGroup group1 = new ButtonGroup();
		for(int i = 0; i < 14; i++){
			JRadioButton button = gp1[i];
			gp1[i].setEnabled(true);
			group1.add(button);
			button.addActionListener(this);
		}
		
		ButtonGroup group2 = new ButtonGroup();
		for(int i = 0; i < 14; i++){
			JRadioButton button = gp2[i];
			gp2[i].setEnabled(false);
			group2.add(button);
			button.addActionListener(this);
		}
		
		ButtonGroup group3 = new ButtonGroup();
		for(int i = 0; i < 14; i++){
			JRadioButton button = gp3[i];
			gp3[i].setEnabled(false);
			group3.add(button);
			button.addActionListener(this);
		}
		
		// Adds actionListeners for the buttons
		JPanel radioPanel = new JPanel(new GridLayout(0,3));
		radioPanel.add(new JLabel("Position 1"));
		radioPanel.add(new JLabel("Position 2"));
		radioPanel.add(new JLabel("Position 3"));
		for(int i = 0; i< 14; i++)
		{
			radioPanel.add(gp1[i]);
			radioPanel.add(gp2[i]);
			radioPanel.add(gp3[i]);
		}
		
		// Adds radioPanel to this RecodeSelector
		add(radioPanel, BorderLayout.CENTER);
	}
	
	/*
		Main action listener
		Invoked when any of the radio buttons is changes
	*/
	public void actionPerformed(ActionEvent e) {
		int selection = Integer.parseInt(e.getActionCommand());
		if(selection <  14){
			current_selection[0] = selection;
		} else if(selection >= 14 && selection < 28){
			current_selection[1] = selection-14;
		} else{
			current_selection[2] = selection-28;
		}
	}
	
	/*
		Returns an array of numbers indicating how the site should be recoded
	*/
	public int[] getRecodeSelection() {
		return current_selection;
	}
	
	/*
		Returns the list of characters of how the sites will be recoded
		Used by the StatsPane.
	*/
	public String getSelectionName(int pos)
	{
		int select = current_selection[pos];
		String name = new String();
		switch(select){
			case 0: name = "A,C,G,T"; break;
			case 1: name = "R,C,T"; break;
			case 2: name = "Y,A,G"; break;
			case 3: name = "R,Y"; break;
			case 4: name = "K,A,C"; break;
			case 5: name = "M,G,T"; break;
			case 6: name = "K,M"; break;
			case 7: name = "S,A,T"; break;
			case 8: name = "W,C,G"; break;
			case 9: name = "S,W"; break;
			case 10: name = "B,A"; break;
			case 11: name = "D,C"; break;
			case 12: name = "H,G"; break;
			case 13: name = "V,T"; break;
		}
		return name;
	}

	/*
		Updates the activity state of the radio buttons to prevent impossibilties
		Called by a BinSelector object when there state is changed.
		For example if the user changes from single nucleotides to dinucleotides
		the second column of recoding (position 2) will become enabled.
	*/
	public void updateState(int bin)
	{
		// If bin selector is set on amino acids - no recoding can occur
		if(bin == 0)
		{
			for(int i = 0; i < 14; i++){
				gp1[i].setEnabled(false);
				gp2[i].setEnabled(false);
				gp3[i].setEnabled(false);
			}
		}
		// If bin selector is set on single nucleotide - only that one position can be recoded
		if(bin == 1)
		{
			for(int i = 0; i < 14; i++){
				gp1[i].setEnabled(true);
				gp2[i].setEnabled(false);
				gp3[i].setEnabled(false);
			}
		}
		// If bin selector is set on duplet nucleotide - positions 1 and 2 can be recoded
		if(bin == 2)
		{
			for(int i = 0; i < 14; i++){
				gp1[i].setEnabled(true);
				gp2[i].setEnabled(true);
				gp3[i].setEnabled(false);
			}
		}
		// If bin selector is set on triplet nucleotide - all positions can be recoded
		if(bin == 3)
		{
			for(int i = 0; i < 14; i++){
				gp1[i].setEnabled(true);
				gp2[i].setEnabled(true);
				gp3[i].setEnabled(true);
			}
		}
	}

}