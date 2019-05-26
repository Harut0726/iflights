package root;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.text.ParseException;
import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.Duration;


public class TimeFormatter {
    public static final float DAY = 1;
    public static final double YEAR = 365.3;

    public static String daysToDate(float days) {
        String oldDate = Params.Mission_Params.START_DATE;
        SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd");
        Calendar c = Calendar.getInstance();
        try {
            c.setTime(sdf.parse(oldDate));
        }
        catch(ParseException e) {
            e.printStackTrace();
        }

        c.add(Calendar.DAY_OF_MONTH, (int)days);
        String newDate = sdf.format(c.getTime());
        return newDate;
    }

    public static long daysBetweenDates(String startDate, String endDate) {
        LocalDate d1 = LocalDate.parse(startDate, DateTimeFormatter.ISO_LOCAL_DATE);
        LocalDate d2 = LocalDate.parse(endDate, DateTimeFormatter.ISO_LOCAL_DATE);
        Duration diff = Duration.between(d1.atStartOfDay(), d2.atStartOfDay());
        long diffDays = diff.toDays();
        return diffDays;
    }
}
