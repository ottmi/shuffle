import java.util.LinkedList;
import java.util.ListIterator;
import javax.swing.*;
import javax.swing.table.*;
import java.awt.GridLayout;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.Component;
import java.awt.*;
import java.lang.Object;
import java.awt.Checkbox;

/*
	Class for a GUI component which allows user selection of individual sequence for
	the purposes of compatibility plot generation
*/
public class SequenceSelector2 extends JPanel{

	private LinkedList<Sequence> seqs;
	private LinkedList<Sequence> returnseqs;
	LinkedList<JCheckBox> check_list;
	SequenceSelectorButtonPanel2 ssbp;
	private int[] selectedRows;
	
	
	/*
		'null' constructor for when no file has been opened
		only appears at program start
	*/
	public SequenceSelector2()
	{
		setLayout(new GridLayout(0,1));
		add(new JLabel("No file selected, please open a file"));
	}
	/*
		Constructs a SequenceSelector for selection of supplied sequences
	*/
	public SequenceSelector2(LinkedList<Sequence> sequences) 
	{
		seqs = sequences;
		setLayout(new BorderLayout());
		check_list= new LinkedList<JCheckBox>();
		ListIterator<Sequence> iter = seqs.listIterator();
		JPanel panel = new JPanel(new GridLayout(0,1));
		while(iter.hasNext())
		{
			String descrip = iter.next().getDescription();
			JCheckBox seqButton = new JCheckBox(descrip);
			seqButton.setSelected(false);
			panel.add(seqButton);
			check_list.add(seqButton);
		}
		
		
		/*
			Listeners for checkbox selection options
		*/
		class SelectAllButtonListener implements ActionListener 
		{
			public void actionPerformed(ActionEvent event) 
			{
				ListIterator<JCheckBox> iter = check_list.listIterator();
				while(iter.hasNext())
				{
					JCheckBox check = iter.next();
					check.setSelected(true);
				}
			}
		}
		class SelectNoneButtonListener implements ActionListener 
		{
			public void actionPerformed(ActionEvent event) 
			{
				ListIterator<JCheckBox> iter = check_list.listIterator();
				while(iter.hasNext())
				{
					JCheckBox check = iter.next();
					check.setSelected(false);
				}
			}
		}
		
		ssbp = new SequenceSelectorButtonPanel2(new SelectAllButtonListener(), new SelectNoneButtonListener());

		//layout and load components

		add(ssbp, BorderLayout.SOUTH);
		JScrollPane panel_t = new JScrollPane(panel);
		add(panel_t, BorderLayout.CENTER);
	}
	
	//returns sequences selected
	public LinkedList<Sequence> getSelectedSequences() 
	{
		LinkedList<Sequence> returnList = new LinkedList<Sequence>();
		int[] selection = getIndices();
		if(selection == null) return new LinkedList<Sequence>();
//		System.out.print("Sequences selected: ");
		for(int i = 0; i < selection.length; i++) {
			Sequence tempseq = seqs.get(selection[i]);
			Sequence newseq = new Sequence(tempseq.getSequence(), tempseq.getDescription());
//			System.out.print(tempseq.getDescription() + ", ");
			returnList.add(newseq);
		}
//		System.out.println("Total = " + returnList.size());
		return returnList;
	}
	
	//helper method, returns indices of selected sequences as
	//they will be ordered the same in the JCheckBox list as in
	//their LinkedList
	private int[] getIndices()
	{
		int size = 0;
		if(check_list == null || check_list.size() == 0) {
			return null;
		}
		ListIterator<JCheckBox> iter = check_list.listIterator();
		while(iter.hasNext())
		{
			JCheckBox check = iter.next();
			if(check.isSelected())
			{
				size++;
			}
		}
		iter = check_list.listIterator();
		int count = 0;
		int i = 0;
		int[] ind = new int[size];
		while(iter.hasNext())
		{
			JCheckBox check = iter.next();
			if(check.isSelected())
			{
				ind[i] = count;
				i++;
			}
			count++;
		}
		return ind;
	}
	
}


/*
	Class for a JPanel containing the buttons for sequence selection options
	Class was made for ease of layout of SequenceSelector component
*/
class SequenceSelectorButtonPanel2 extends JPanel
{
	static final long serialVersionUID = 1125733L;
	
	public SequenceSelectorButtonPanel2(ActionListener sal, ActionListener snl) {

		JButton selectAllButton = new JButton("Select All");
		JButton selectNoneButton = new JButton("Select None");
		

		selectAllButton.addActionListener(sal);
		selectNoneButton.addActionListener(snl);
		
		setLayout(new GridLayout(1,2));
		add(selectAllButton);
		add(selectNoneButton);
	}
}