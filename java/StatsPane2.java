import javax.swing.*;

/*
	This class describes an object which is essentially a JTextField inside a JScrollPane, which is used
	in lieu of the console for the purposes of feedback of information to the user. It automatically scrolls
	to the bottom of the pane when text is added.
*/
public class StatsPane2 extends JScrollPane {
	JTextArea area;
	
	public StatsPane2() {
		super();
		area = new JTextArea();
		setViewportView(area);
		area.setLineWrap(true);
		area.setEditable(false);
	}
	
	public void append(String s) {
//		area.append(s);
		
		// Determine whether the scrollbar is currently at the very bottom position.
		JScrollBar vbar = getVerticalScrollBar();
		boolean autoScroll = ((vbar.getValue() + vbar.getVisibleAmount()) == vbar.getMaximum());
		
		// append to the JTextArea (that's wrapped in a JScrollPane named 'scrollPane'
		area.append(s);

		// now scroll if we were already at the bottom.
		if( autoScroll ) area.setCaretPosition( area.getDocument().getLength() );
	}
	
	public void appendln(String s) {
		append(s + "\n");		
	}
	
}