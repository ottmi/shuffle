import java.awt.*;
import javax.swing.*;
import java.awt.event.*;

public class BinSelector extends JPanel implements ActionListener {
	
	private int bin_selection_number;
	private int bin_offset_number;
	private int grouping_number;
	private JPanel radioPanel;
	private JRadioButton aminoAcids, binOne, binTwo, binThree, standard, offsetZero, offsetOne, offsetTwo, groupYes, groupNo;
	private RecodeSelector recode;
	private SiteSelector select;
	private boolean group;
	
	public BinSelector(RecodeSelector rs, SiteSelector ss) {
		super(new BorderLayout());
		
		bin_selection_number = 1;	//default value
		bin_offset_number = 0;	//default value
		grouping_number = 1;
		recode = rs;
		select=ss;
		group = false;
		
		//creates radio buttons for bin grouping number
		
		aminoAcids = new JRadioButton("Amino Acids");
		aminoAcids.setActionCommand("0");
		
		binOne = new JRadioButton("Nucleotide (singlets)");
		binOne.setActionCommand("1");
		binOne.setSelected(true);
		
		binTwo = new JRadioButton("Nucleotide (duplets)");
		binTwo.setActionCommand("2");
	
		
		binThree = new JRadioButton("Nucleotide (triplets)");
		binThree.setActionCommand("3");
		
		standard = new JRadioButton("Standard (A-Z, 0-9)");
		standard.setActionCommand("4");
		
		
		
		//creates radio buttons for bin offset number
		
		offsetZero = new JRadioButton("Frame 1");
		offsetZero.setActionCommand("0");
		offsetZero.setSelected(true);
		offsetZero.setEnabled(false);
		
		offsetOne = new JRadioButton("Frame 2");
		offsetOne.setActionCommand("1");
		offsetOne.setEnabled(false);
		
		offsetTwo = new JRadioButton("Frame 3");
		offsetTwo.setActionCommand("2");
		offsetTwo.setEnabled(false);
		
		//creates radio buttons for grouping
		
		groupYes = new JRadioButton("Yes");
		groupYes.setActionCommand("0");
		groupYes.setSelected(false);
		groupYes.setEnabled(true);
		
		groupNo = new JRadioButton("No");
		groupNo.setActionCommand("1");
		groupNo.setSelected(true);
		groupNo.setEnabled(true);
		
		
		//button group for bin grouping selectors
		ButtonGroup group = new ButtonGroup();
		group.add(aminoAcids);
		group.add(binOne);
		group.add(binTwo);
		group.add(binThree);
		group.add(standard);
		
		//button group for bin offset selectors
		ButtonGroup group_offset = new ButtonGroup();
		group_offset.add(offsetZero);
		group_offset.add(offsetOne);
		group_offset.add(offsetTwo);
		
		//button group for grouping
		ButtonGroup group_grouping = new ButtonGroup();
		group_grouping.add(groupYes);
		group_grouping.add(groupNo);
		
		
		//set ActionListeners for all of the bin grouping radiobuttons
		aminoAcids.addActionListener(this);
		binOne.addActionListener(this);
		binTwo.addActionListener(this);
		binThree.addActionListener(this);
		standard.addActionListener(this);
		
		//define an actionlistener for offset buttons, add listener to offset buttons
		class BinOffsetListener implements ActionListener {
			public void actionPerformed(ActionEvent e) {
				bin_offset_number = Integer.parseInt(e.getActionCommand());
			}
		}
		
		offsetZero.addActionListener(new BinOffsetListener());
		offsetOne.addActionListener(new BinOffsetListener());
		offsetTwo.addActionListener(new BinOffsetListener());
		
		//define an actionlistener for grouping buttons, add listener to grouping buttons
		class GroupingListener implements ActionListener {
			public void actionPerformed(ActionEvent e) {
				grouping_number = Integer.parseInt(e.getActionCommand());
				select.updateState(bin_selection_number,grouping_number);
			}
		}
		
		groupYes.addActionListener(new GroupingListener());
		groupNo.addActionListener(new GroupingListener());
		
		//puts the radio buttons in a column in a panel
		radioPanel = new JPanel(new GridLayout(15,1));
		radioPanel.add(new JLabel("Data Type"));
		radioPanel.add(aminoAcids);
		radioPanel.add(binOne);
		radioPanel.add(binTwo);
		radioPanel.add(binThree);
		radioPanel.add(standard);
		radioPanel.add(new JLabel("Frames"));
		radioPanel.add(offsetZero);
		radioPanel.add(offsetOne);
		radioPanel.add(offsetTwo);
		radioPanel.add(new JLabel("Group?"));
		radioPanel.add(groupYes);
		radioPanel.add(groupNo);
		
		//adds radioPanel to the BinSelector
		add(radioPanel, BorderLayout.CENTER);
	}
	
	public void updateState()
	{
		if(bin_selection_number == 0)
		{
			offsetZero.setEnabled(false);
			offsetOne.setEnabled(false);
			offsetTwo.setEnabled(false);
			groupYes.setEnabled(false);
			groupNo.setEnabled(false);
		}
		if(bin_selection_number == 1)
		{
			offsetZero.setEnabled(false);
			offsetOne.setEnabled(false);
			offsetTwo.setEnabled(false);
			groupYes.setEnabled(false);
			groupNo.setEnabled(false);
		}
		if(bin_selection_number == 2)
		{
			offsetZero.setEnabled(true);
			offsetOne.setEnabled(true);
			offsetTwo.setEnabled(false);
			groupYes.setEnabled(true);
			groupNo.setEnabled(true);
		}
		if(bin_selection_number == 3)
		{
			offsetZero.setEnabled(true);
			offsetOne.setEnabled(true);
			offsetTwo.setEnabled(true);		
			groupYes.setEnabled(true);
			groupNo.setEnabled(true);
		}
		if(bin_selection_number == 4)
		{
			offsetZero.setEnabled(false);
			offsetOne.setEnabled(false);
			offsetTwo.setEnabled(false);
			groupYes.setEnabled(false);
			groupNo.setEnabled(false);		
		}
	}
	
	public void actionPerformed(ActionEvent e)
	{
		bin_selection_number = Integer.parseInt(e.getActionCommand());
		updateState();
		recode.updateState(bin_selection_number);
		select.updateState(bin_selection_number,grouping_number);
	}
	
	public int getBinSelectionNumber() {return bin_selection_number; }
	
	public int getBinOffsetNumber() { return bin_offset_number; }
	
	public int getGrouping() { return grouping_number; }

}