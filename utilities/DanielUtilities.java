package utilities;

import java.util.*;

/**
 * This has a lot of utility methods that I'll use to supplement the rest of my code. Some of this
 * code will be taken from the Berkeley parser and other sources (and I'll cite them, of course).
 * 
 * @author Daniel Seita
 */
public class DanielUtilities {
    
    /**
     * Given log(x) and log(y), returns log(x+y). This method is from Dan Klein, but the formula is
     * pretty commonplace.
     * 
     * @param logX This is log(x) for some x.
     * @param logY This is log(y) for some y.
     * @return The value log(x+y).
     */
	public static double logAdd(double logX, double logY) {
		if (logY > logX) {
			double temp = logX;
			logX = logY;
			logY = temp;
		}
		if (logX == Double.NEGATIVE_INFINITY) { return logX; }
		double negDiff = logY - logX;
		if (negDiff < -20) { return logX; }
		return logX + java.lang.Math.log(1.0 + java.lang.Math.exp(negDiff));
	}

	/**
	 * Simple method which turns an array of command line arguments into a map, where each token
	 * starting with a '-' is a key and the following non '-' initial token, if there is one, is
	 * the value. For example, '-size 5 -verbose' will produce keys (-size,5) and (-verbose,null).
	 * This method was written by Dan Klein.
	 */
	public static Map<String, String> simpleCommandLineParser(String[] args) {
		Map<String, String> map = new HashMap<String, String>();
		for (int i = 0; i <= args.length; i++) {
			String key = (i > 0 ? args[i - 1] : null);
			String value = (i < args.length ? args[i] : null);
			if (key == null || key.startsWith("-")) {
				if (value != null && value.startsWith("-")) value = null;
				if (key != null || value != null) map.put(key, value);
			}
		}
		return map;
	}

	/**
	 * Simple method to look up a key in an argument map. Returns the defaultValue if the argument
	 * is not specified in the map. This method was written by Dan Klein.
	 */
	public static String getValueOrUseDefault(Map<String, String> args, String key, String defaultValue) {
		if (args.containsKey(key)) return args.get(key);
		return defaultValue;
	}
}