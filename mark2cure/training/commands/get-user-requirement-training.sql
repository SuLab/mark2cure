SELECT  `task_requirement`.`hash`,
        ANY_VALUE(`task_requirement`.`name`) as `name`,
        MAX(`task_level`.`created`) as `last_created`,
        COUNT(`task_level`.`id`) as `completions`

FROM `task_requirement`

LEFT JOIN `task_level`
  ON `task_level`.`requirement_id` = `task_requirement`.`id`
  AND `task_level`.`user_id` = {user_id}

WHERE `task_requirement`.`task_type` = "{task_type}"

GROUP BY `task_requirement`.`hash`, `task_requirement`.`order`

ORDER BY `task_requirement`.`order` ASC
