import java.awt.*;
import javax.swing.*;
import java.lang.*;

/*
	Trivial class for display of application
*/
public class Shuffle {
	
	public static void main(String[] args) {
		JFrame appFrame = new ApplicationGui("Shuffle");
		appFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		appFrame.pack();
		appFrame.setVisible(true);
	}
}