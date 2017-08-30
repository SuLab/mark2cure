SELECT  `task_level`.`task_type`,
        `task_level`.`level`,
        MAX(`task_level`.`created`) as `last_created`,
        COUNT(`task_level`.`id`) as `completions`
FROM `task_level`
WHERE `task_level`.`user_id` = {user_id}
GROUP BY `task_type` ASC, `level`
ORDER BY `task_type` ASC, `level`
