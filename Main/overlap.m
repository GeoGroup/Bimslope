function result = overlap(a, b)
    if mycontains(a(1), b)
        result = true;
        return
    end
    if mycontains(a(2), b)
        result = true;
        return
    end
    if mycontains(b(1), a)
        result = true;
        return
    end
    if mycontains(b(2), a)
        result = true;
        return
    end

    result = false;
end


