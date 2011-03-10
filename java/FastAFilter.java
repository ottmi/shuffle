import java.io.File;
import javax.swing.*;
import javax.swing.filechooser.*;

/*
	Filter for JFileChooser
*/
public class FastAFilter extends FileFilter {
	public boolean accept(File f) {
		return true;
	}
	public String getDescription() {
		return "FastA Files";
	}
}